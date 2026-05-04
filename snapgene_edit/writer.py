"""Binary writer for SnapGene .dna files.

Preserves the original block order for maximum compatibility with SnapGene.
"""

from __future__ import annotations

import struct
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

from snapgene_edit.models import Feature, SnapGeneFile


def write_snapgene(sgf: SnapGeneFile, path: str) -> None:
    """Write a SnapGeneFile to a .dna file.

    Args:
        sgf: The SnapGeneFile object to write.
        path: Path to the output .dna file.
    """
    buffer = _serialize(sgf)
    with open(path, "wb") as f:
        f.write(buffer)


def _make_block(block_type: int, data: bytes) -> bytes:
    """Build a binary block: type(1 byte) + size(4 bytes big-endian) + data(N bytes)."""
    return struct.pack("B", block_type) + struct.pack(">I", len(data)) + data


def _make_seq_block(sgf: SnapGeneFile) -> bytes:
    """Build the DNA sequence block (type 0)."""
    seq_data = sgf.sequence.encode("ascii")
    flags = 0
    if sgf.is_circular:
        flags |= 0x01
    if sgf.is_double_stranded:
        flags |= 0x02
    if sgf.dam_methylated:
        flags |= 0x04
    if sgf.dcm_methylated:
        flags |= 0x08
    if sgf.eco_ki_methylated:
        flags |= 0x10
    return _make_block(0, bytes([flags]) + seq_data)


def _serialize(sgf: SnapGeneFile) -> bytes:
    """Serialize a SnapGeneFile to binary format, preserving original block order."""
    # Build header
    header = b"".join([
        b"\t",                         # first byte (tab)
        struct.pack(">I", 14),         # length
        b"SnapGene",                   # title
        struct.pack(">H", 1 if sgf.is_dna else 0),
        struct.pack(">H", sgf.export_version),
        struct.pack(">H", sgf.import_version),
    ])

    # Build replacement blocks for types we know how to regenerate
    replacements: Dict[int, bytes] = {}

    # Block 0: DNA Sequence
    replacements[0] = _make_seq_block(sgf)

    # Block 6: Notes - only replace if we have notes that differ from raw
    if sgf.notes:
        raw_notes = _get_raw_block(sgf._raw_blocks, 6)
        new_notes_xml = _generate_notes_xml(sgf.notes).encode("utf-8")
        # If the raw notes match, keep them to avoid changing binary output
        if raw_notes != new_notes_xml:
            replacements[6] = _make_block(6, new_notes_xml)
        elif raw_notes is not None:
            replacements[6] = _make_block(6, raw_notes)
    else:
        raw_notes = _get_raw_block(sgf._raw_blocks, 6)
        if raw_notes is not None:
            replacements[6] = _make_block(6, raw_notes)

    # Block 10: Features - only replace if we have features
    if sgf.features:
        raw_feats = _get_raw_block(sgf._raw_blocks, 10)
        new_feats_xml = _generate_features_xml(sgf.features).encode("utf-8")
        if raw_feats != new_feats_xml:
            replacements[10] = _make_block(10, new_feats_xml)
        elif raw_feats is not None:
            replacements[10] = _make_block(10, raw_feats)

    # Block 5: Primers - keep raw unless we have primer edits
    raw_primer = _get_raw_block(sgf._raw_blocks, 5)
    if raw_primer is not None:
        if sgf.primers:
            # We have primers data but no primer XML generation yet
            # For now, just use raw data
            pass
        replacements[5] = _make_block(5, raw_primer)

    # Iterate original raw blocks IN ORDER and replace where needed
    parts: List[bytes] = [header]
    written_types: set = set()

    for raw in sgf._raw_blocks:
        bt = raw["type"]
        if bt in replacements:
            parts.append(replacements[bt])
            written_types.add(bt)
        else:
            parts.append(_make_block(bt, raw["data"]))
            written_types.add(bt)

    # Write any replacement blocks that weren't in the original order
    for bt in sorted(replacements.keys()):
        if bt not in written_types:
            parts.append(replacements[bt])
            written_types.add(bt)

    return b"".join(parts)


def _generate_features_xml(features: List[Feature]) -> str:
    """Generate the Features XML string from Feature objects."""
    root = ET.Element("Features")

    for feat in features:
        feat_elem = ET.SubElement(root, "Feature")
        feat_elem.set("name", feat.name)
        feat_elem.set("type", feat.type)
        dir_map = {0: "0", 1: "1", -1: "2"}
        feat_elem.set("directionality", dir_map.get(feat.strand, "0"))

        if feat.segments:
            for seg in feat.segments:
                seg_elem = ET.SubElement(feat_elem, "Segment")
                seg_elem.set("range", f"{seg.range_start}-{seg.range_end}")
                seg_elem.set("color", seg.color)
                seg_elem.set("type", seg.segment_type)
                if seg.name:
                    seg_elem.set("name", seg.name)
        else:
            seg_elem = ET.SubElement(feat_elem, "Segment")
            seg_elem.set("range", f"{feat.start + 1}-{feat.end}")
            seg_elem.set("color", feat.color)
            seg_elem.set("type", "standard")

        for q_name, q_value in feat.qualifiers.items():
            q_elem = ET.SubElement(feat_elem, "Q")
            q_elem.set("name", q_name)

            values = q_value if isinstance(q_value, list) else [q_value]
            for v in values:
                v_elem = ET.SubElement(q_elem, "V")
                if isinstance(v, int):
                    v_elem.set("int", str(v))
                    v_elem.text = str(v)
                elif isinstance(v, float):
                    v_elem.set("float", str(v))
                    v_elem.text = str(v)
                else:
                    v_elem.text = str(v)

    return ET.tostring(root, encoding="unicode")


def _generate_notes_xml(notes: Dict[str, Any]) -> str:
    """Generate the Notes XML string from a dictionary.

    Supports dicts with ``_text`` key for elements that have both
    text content and attributes (e.g., ``<Created UTC="...">date</Created>``).
    """
    root = ET.Element("Notes")

    for key, value in notes.items():
        child = ET.SubElement(root, key)
        if isinstance(value, dict) and "_text" in value:
            # Element with both text content and attributes
            child.text = str(value["_text"])
            for attr_key, attr_val in value.items():
                if attr_key != "_text":
                    child.set(attr_key, str(attr_val))
        elif isinstance(value, dict):
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, dict) and "_text" in sub_value:
                    # Sub-element with text + attributes
                    sub_elem = ET.SubElement(child, sub_key)
                    sub_elem.text = str(sub_value["_text"])
                    for attr_key, attr_val in sub_value.items():
                        if attr_key != "_text":
                            sub_elem.set(attr_key, str(attr_val))
                else:
                    sub_elem = ET.SubElement(child, sub_key)
                    sub_elem.text = str(sub_value)
        elif isinstance(value, list):
            for item in value:
                item_elem = ET.SubElement(child, "Item")
                item_elem.text = str(item)
        else:
            child.text = str(value)

    return ET.tostring(root, encoding="unicode")


def _get_raw_block(raw_blocks: List[Dict[str, Any]], block_type: int) -> Optional[bytes]:
    """Get raw data for a specific block type from preserved blocks."""
    for raw in raw_blocks:
        if raw["type"] == block_type:
            return raw["data"]
    return None
