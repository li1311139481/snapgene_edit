"""Binary reader for SnapGene .dna files."""

from __future__ import annotations

import struct
import xml.etree.ElementTree as ET
from typing import Any, Dict, List

from snapgene_edit.models import Feature, FeatureSegment, SnapGeneFile


def read_snapgene(path: str) -> SnapGeneFile:
    """Read a SnapGene .dna file and return a SnapGeneFile object.

    Args:
        path: Path to the .dna file.

    Returns:
        A SnapGeneFile object with parsed contents.
    """
    with open(path, "rb") as f:
        data = f.read()
    return _parse_buffer(data, path)


def _parse_buffer(buffer: bytes, source_name: str = "") -> SnapGeneFile:
    """Parse a binary buffer containing SnapGene format data."""
    offset = 0

    def read(size: int) -> bytes:
        nonlocal offset
        start = offset
        offset += size
        return buffer[start:offset]

    # --- Parse header ---
    read(1)  # first byte (tab character 0x09), skip it

    length = struct.unpack(">I", read(4))[0]
    title = read(8).decode("ascii")

    if length != 14 or title != "SnapGene":
        raise ValueError(
            f"Wrong format for a SnapGene file: length={length}, title='{title}'"
        )

    is_dna = struct.unpack(">H", read(2))[0]
    export_version = struct.unpack(">H", read(2))[0]
    import_version = struct.unpack(">H", read(2))[0]

    result = SnapGeneFile(
        name=source_name,
        is_dna=bool(is_dna),
        export_version=export_version,
        import_version=import_version,
    )

    raw_blocks: List[Dict[str, Any]] = []

    # --- Parse blocks ---
    while offset < len(buffer):
        next_byte = read(1)
        block_size = struct.unpack(">I", read(4))[0]
        block_type = next_byte[0]

        # Always store raw block data for preservation
        raw_data = read(block_size)

        if block_type == 0:
            # DNA Sequence block
            flags = raw_data[0]
            seq_len = block_size - 1
            if seq_len < 0:
                raise ValueError("Failed parsing SnapGene: < 0 length sequence")
            seq_bytes = raw_data[1:]
            result.sequence = seq_bytes.decode("ascii")
            result.is_circular = bool(flags & 0x01)
            result.is_double_stranded = bool(flags & 0x02)
            result.dam_methylated = bool(flags & 0x04)
            result.dcm_methylated = bool(flags & 0x08)
            result.eco_ki_methylated = bool(flags & 0x10)

        elif block_type == 5:
            # Primers
            pass

        elif block_type == 6:
            # Notes (XML)
            xml_data = raw_data.decode("utf-8")
            result.notes = _parse_notes_xml(xml_data)

        elif block_type == 8:
            # Additional sequence properties
            pass

        elif block_type == 10:
            # Features (XML)
            xml_data = raw_data.decode("utf-8")
            result.features = _parse_features_xml(xml_data)

        elif block_type == 11:
            # History node
            pass

        # Always store the raw block for re-serialization order preservation
        raw_blocks.append({
            "type": block_type,
            "data": raw_data,
        })

    result._raw_blocks = raw_blocks
    return result


def _parse_features_xml(xml_str: str) -> List[Feature]:
    """Parse the Features XML block into Feature objects."""
    features: List[Feature] = []
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError:
        return features

    for feat_elem in root.findall("Feature"):
        name = feat_elem.get("name", "Unnamed Feature")
        ftype = feat_elem.get("type", "misc_feature")
        dir_str = feat_elem.get("directionality", "0")
        strand_map = {"0": 0, "1": 1, "2": -1, "3": 0}
        strand = strand_map.get(dir_str, 0)

        # Parse segments
        segment_elems = feat_elem.findall("Segment")
        segments: List[FeatureSegment] = []
        min_start = 0
        max_end = 0

        for seg_elem in segment_elems:
            range_str = seg_elem.get("range", "1-1")
            parts = range_str.split("-")
            seg_start = int(parts[0])
            seg_end = int(parts[1])
            seg_color = seg_elem.get("color", "#999999")
            seg_type = seg_elem.get("type", "standard")
            seg_name = seg_elem.get("name", "")

            if min_start == 0 or seg_start < min_start:
                min_start = seg_start
            if seg_end > max_end:
                max_end = seg_end

            segments.append(FeatureSegment(
                range_start=seg_start,
                range_end=seg_end,
                color=seg_color,
                segment_type=seg_type,
                name=seg_name,
            ))

        if min_start == 0 and max_end == 0:
            min_start = 1
            max_end = 1

        # Parse qualifiers
        qualifiers: Dict[str, Any] = {}
        for q_elem in feat_elem.findall("Q"):
            q_name = q_elem.get("name", "")
            v_elems = q_elem.findall("V")
            values = []
            for v in v_elems:
                if v.text:
                    values.append(v.text)
                # Check for attributes on V element
                for attr_key in ("text", "int", "float", "bool"):
                    attr_val = v.get(attr_key)
                    if attr_val is not None:
                        if attr_key == "int":
                            values.append(int(attr_val))
                        elif attr_key == "float":
                            values.append(float(attr_val))
                        elif attr_key == "bool":
                            values.append(attr_val.lower() == "true")
                        else:
                            values.append(attr_val)
            if len(values) == 1:
                qualifiers[q_name] = values[0]
            elif len(values) > 1:
                qualifiers[q_name] = values

        # Use first segment's color if available
        color = segments[0].color if segments else "#999999"

        feature = Feature(
            name=name,
            type=ftype,
            start=min_start - 1,
            end=max_end,
            strand=strand,
            color=color,
            segments=segments,
            qualifiers=qualifiers,
        )
        features.append(feature)

    return features


def _parse_notes_xml(xml_str: str) -> Dict[str, Any]:
    """Parse the Notes XML block into a dictionary.

    Elements with attributes are stored as dicts with a ``_text`` key
    for the text content and other keys for attributes.
    """
    notes: Dict[str, Any] = {}
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError:
        return notes

    for child in root:
        tag = child.tag
        text = child.text.strip() if child.text else ""
        attribs = dict(child.attrib)
        sub_elems = list(child)

        if attribs:
            # Element has attributes → store as dict
            entry: Dict[str, Any] = {}
            if text:
                entry["_text"] = text
            for k, v in attribs.items():
                entry[k] = v
            if sub_elems:
                for sub in sub_elems:
                    sub_text = sub.text.strip() if sub.text else ""
                    entry[sub.tag] = sub_text
            notes[tag] = entry

        elif sub_elems:
            # Element has sub-elements (but no attributes on itself)
            sub_dict = {}
            for sub in sub_elems:
                sub_tag = sub.tag
                sub_text = sub.text.strip() if sub.text else ""
                if sub.attrib:
                    attr_entry: Dict[str, Any] = {}
                    if sub_text:
                        attr_entry["_text"] = sub_text
                    for k, v in sub.attrib.items():
                        attr_entry[k] = v
                    sub_dict[sub_tag] = attr_entry
                else:
                    sub_dict[sub_tag] = sub_text
            notes[tag] = sub_dict

        elif text:
            notes[tag] = text
        else:
            notes[tag] = ""

    return notes
