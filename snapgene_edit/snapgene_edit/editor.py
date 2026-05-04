"""High-level editor API for SnapGene files."""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from snapgene_edit.reader import read_snapgene
from snapgene_edit.writer import write_snapgene
from snapgene_edit.models import SnapGeneFile, Feature, FeatureSegment
from snapgene_edit.config import FLUORESCENT_FILES
from snapgene_edit.utils import reverse_complement as revcomp


class SnapGeneEditor:
    """High-level editor for SnapGene files.

    Provides methods to load, edit, and save SnapGene .dna files.

    Usage:
        editor = SnapGeneEditor("plasmid.dna")
        print(editor.sequence)
        editor.rename_feature("old_name", "new_name")
        editor.add_feature("GFP", 100, 800, type="CDS")
        editor.save("modified.dna")
    """

    def __init__(self, filepath: Optional[str] = None, data: Optional[SnapGeneFile] = None):
        """Initialize the editor.

        Args:
            filepath: Path to a .dna file to load.
            data: An existing SnapGeneFile object.
        """
        if filepath:
            self._filepath = filepath
            self._data = read_snapgene(filepath)
        elif data:
            self._data = data
            self._filepath = ""
        else:
            self._data = SnapGeneFile()
            self._filepath = ""

    # --- Convenience properties ---

    @property
    def name(self) -> str:
        return self._data.name

    @name.setter
    def name(self, value: str) -> str:
        self._data.name = value

    @property
    def sequence(self) -> str:
        return self._data.sequence

    @property
    def length(self) -> int:
        return len(self._data.sequence)

    @property
    def features(self) -> List[Feature]:
        return self._data.features

    @property
    def is_circular(self) -> bool:
        return self._data.is_circular

    @is_circular.setter
    def is_circular(self, value: bool) -> str:
        self._data.is_circular = value

    # --- Feature coordinate adjustment ---

    def _adjust_features_for_edit(self, edit_start: int, edit_end: int,
                                   old_len: int, new_seq: str) -> None:
        """Adjust feature coordinates after a sequence edit.

        Rules:
          - Feature entirely **before** the edit region → unchanged.
          - Feature entirely **after** the edit region → both ``start`` and
            ``end`` are shifted by ``len(new_seq) - old_len``.
          - Feature **completely within** ``[edit_start, edit_end)`` →
            sequence fully changed → feature is deleted.
          - Feature **partially overlapping** the edit region →
            the overlapping portion is deleted; the non-overlapping portion(s)
            survive with their **original name**.

        Args:
            edit_start: 0-based start of the edit.
            edit_end: End of the edit *before* the edit (exclusive).
            old_len: Original length of the replaced region (edit_end - edit_start).
            new_seq: The new sequence string being inserted.
        """
        delta = len(new_seq) - old_len
        new_features = []

        for f in self._data.features:
            fs, fe = f.start, f.end

            # --- Case 1: entirely before the edit → unchanged ---
            if fe <= edit_start:
                new_features.append(f)
                continue

            # --- Case 2: entirely after the edit → shift by delta ---
            if fs >= edit_end:
                f.start = fs + delta
                f.end = fe + delta
                for seg in f.segments:
                    seg.range_start += delta
                    seg.range_end += delta
                new_features.append(f)
                continue

            # --- Case 3 & 4: overlaps the edit region ---
            # Determine which non-overlapping parts survive
            has_left = fs < edit_start      # [fs, edit_start) survives
            has_right = fe > edit_end        # [edit_end, fe) → shifted, survives

            if not has_left and not has_right:
                # Case 3: completely contained → delete (don't add)
                continue

            # Case 4: partially overlapping
            # Each surviving part keeps the original name (no _partion suffix)

            if has_left and not has_right:
                # Only left part survives → trim end to edit_start
                f.end = edit_start
                if f.segments:
                    f.segments[0].range_end = edit_start
                new_features.append(f)

            elif has_right and not has_left:
                # Only right part survives → shift start/end after the new sequence
                new_start = edit_start + len(new_seq)
                orig_len = fe - edit_end
                f.start = new_start
                f.end = new_start + orig_len
                if f.segments:
                    f.segments[0].range_start = new_start + 1
                    f.segments[0].range_end = new_start + orig_len
                new_features.append(f)

            else:
                # has_left and has_right: both survive
                # Keep original feature for the left part, truncated
                f.end = edit_start
                if f.segments:
                    f.segments[0].range_end = edit_start
                new_features.append(f)

                # Create a new feature for the right part (same name)
                new_start = edit_start + len(new_seq)
                orig_len = fe - edit_end
                right = Feature(
                    name=f.name,
                    type=f.type,
                    start=new_start,
                    end=new_start + orig_len,
                    strand=f.strand,
                    color=f.color,
                    qualifiers=f.qualifiers.copy(),
                    segments=[FeatureSegment(
                        range_start=new_start + 1,
                        range_end=new_start + orig_len,
                        color=f.color,
                    )],
                )
                new_features.append(right)

        self._data.features = new_features

    # --- Sequence editing by position ---

    def set_sequence(self, sequence: str) -> str:
        """Replace the entire sequence (clears all feature annotations)."""
        self._data.sequence = sequence.upper()
        self._data.features = []

    def replace_region(self, start: int, end: int, new_seq: str) -> str:
        """Replace a region of the sequence (0-based, end is exclusive).

        Features are automatically adjusted:
          - Features before the replaced region → unchanged.
          - Features after the replaced region → shifted by the length delta.
          - Features completely within the replaced region → deleted.
          - Features partially overlapping → the overlapping portion is
            deleted; the non-overlapping portion survives with its original name.
        """
        seq = self._data.sequence
        if start < 0 or end > len(seq) or start > end:
            raise ValueError(
                f"Invalid range: ({start}, {end}) for sequence length {len(seq)}"
            )
        old_len = end - start
        new_seq_upper = new_seq.upper()
        self._adjust_features_for_edit(start, end, old_len, new_seq_upper)
        self._data.sequence = seq[:start] + new_seq_upper + seq[end:]

    def delete_region(self, start: int, end: int) -> str:
        """Delete a region of the sequence (0-based, end is exclusive).

        Features are automatically adjusted:
          - Features before the deleted region → unchanged.
          - Features after the deleted region → shifted left.
          - Features completely within the deleted region → deleted.
          - Features partially overlapping → the overlapping portion is
            deleted; the non-overlapping portion survives with its original name.
        """
        self.replace_region(start, end, "")

    def insert(self, position: int, new_seq: str) -> str:
        """Insert a sequence at a given position (0-based).

        Features are automatically adjusted:
          - Features before the insertion point → unchanged.
          - Features at or after the insertion point → shifted right.
          - Features straddling the insertion point → the overlapping portion
            is deleted; the non-overlapping portion survives with its
            original name.
        """
        seq = self._data.sequence
        if position < 0 or position > len(seq):
            raise ValueError(
                f"Invalid position: {position} for sequence length {len(seq)}"
            )
        new_seq_upper = new_seq.upper()
        self._adjust_features_for_edit(position, position, 0, new_seq_upper)
        self._data.sequence = seq[:position] + new_seq_upper + seq[position:]

    # --- Sequence editing by feature name or sequence search ---

    def replace_in_feature(self, feature_name: str, new_seq: str) -> bool:
        """Replace the entire region of a named feature with new_seq.

        The target feature itself is preserved (its coordinate is adjusted
        to match the new sequence length). All other overlapping features
        are split into ``{name}_partion`` parts.

        Returns True if the feature was found and replaced.
        """
        feat = self._data.find_feature(feature_name)
        if not feat:
            return False

        old_start, old_end = feat.start, feat.end
        new_seq_upper = new_seq.upper()

        # Temporarily remove the target feature so it's not split by replace_region
        self._data.remove_feature(feature_name)

        # Replace the region (handles all OTHER features)
        self.replace_region(old_start, old_end, new_seq_upper)

        # Re-add the target feature with adjusted coordinates
        feat.start = old_start
        feat.end = old_start + len(new_seq_upper)
        if feat.segments:
            feat.segments[0].range_start = feat.start + 1
            feat.segments[0].range_end = feat.end
        self._data.add_feature(feat)
        return True

    def insert_in_feature(self, feature_name: str, new_seq: str,
                          *, at: str = "start") -> bool:
        """Insert new_seq inside a named feature.

        The target feature itself is preserved (its end coordinate is
        extended). All other overlapping features are split into
        ``{name}_partion`` parts.

        Args:
            feature_name: Name of the feature.
            new_seq: Sequence to insert.
            at: Where to insert - "start" (default), "end", or "center".

        Returns True if the feature was found and the insertion was done.
        """
        feat = self._data.find_feature(feature_name)
        if not feat:
            return False

        if at == "start":
            pos = feat.start
        elif at == "end":
            pos = feat.end
        elif at == "center":
            pos = (feat.start + feat.end) // 2
        else:
            raise ValueError(
                f"Invalid 'at' value: {at!r}. Use 'start', 'end', or 'center'."
            )

        # Temporarily remove the target feature so it's not split by insert
        self._data.remove_feature(feature_name)

        # Insert at position (handles all OTHER features)
        self.insert(pos, new_seq)

        # Re-add the target feature with extended coordinates
        new_seq_upper = new_seq.upper()
        feat.end += len(new_seq_upper)
        if feat.segments:
            feat.segments[0].range_start = feat.start + 1
            feat.segments[0].range_end = feat.end
        self._data.add_feature(feat)
        return True

    def find_sequence(self, seq: str, *, case_sensitive: bool = False) -> List[int]:
        """Find all occurrences of seq in the current sequence.

        Args:
            seq: Subsequence to search for.
            case_sensitive: Whether the search is case-sensitive.

        Returns:
            List of 0-based start positions where seq is found.
        """
        s = self._data.sequence
        if not case_sensitive:
            s = s.upper()
            search = seq.upper()
        else:
            search = seq
        positions = []
        start = 0
        while True:
            pos = s.find(search, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions

    def replace_sequence(self, old_seq: str, new_seq: str, *, all: bool = True) -> int:
        """Replace occurrences of old_seq with new_seq in the current sequence.

        Args:
            old_seq: Sequence to find and replace.
            new_seq: Replacement sequence.
            all: If True (default), replace all occurrences.
                 If False, replace only the first occurrence.

        Returns:
            Number of occurrences replaced.
        """
        positions = self.find_sequence(old_seq, case_sensitive=False)
        if not positions:
            return 0

        if all:
            positions = sorted(positions, reverse=True)
            for pos in positions:
                self.replace_region(pos, pos + len(old_seq), new_seq)
        else:
            pos = positions[0]
            self.replace_region(pos, pos + len(old_seq), new_seq)
        return len(positions)

    # --- Feature editing ---

    def add_feature(
        self,
        name: str,
        start: int,
        end: int,
        *,
        feature_type: str = "misc_feature",
        strand: int = 0,
        color: str = "#999999",
    ) -> Feature:
        """Add a new feature.

        Args:
            name: Feature name.
            start: 0-based start position.
            end: 0-based end position (exclusive).
            feature_type: Feature type (e.g., "CDS", "gene", "promoter").
            strand: 1 (forward), -1 (reverse), 0 (none).
            color: Hex color string.

        Returns:
            The created Feature.
        """
        feature = Feature(
            name=name,
            type=feature_type,
            start=start,
            end=end,
            strand=strand,
            color=color,
            segments=[FeatureSegment(range_start=start + 1, range_end=end, color=color)],
        )
        self._data.add_feature(feature)
        return feature

    def remove_feature(self, name: str) -> bool:
        """Remove a feature by name.

        Returns True if the feature was found and removed.
        """
        before = len(self._data.features)
        self._data.remove_feature(name)
        return len(self._data.features) < before

    def rename_feature(self, old_name: str, new_name: str) -> bool:
        """Rename a feature.

        Returns True if the feature was found and renamed.
        """
        feat = self._data.find_feature(old_name)
        if feat:
            feat.name = new_name
            return True
        return False

    def update_feature_coords(self, name: str, start: int, end: int) -> bool:
        """Update the coordinates of a feature.

        Returns True if the feature was found and updated.
        """
        feat = self._data.find_feature(name)
        if feat:
            feat.start = start
            feat.end = end
            if feat.segments:
                feat.segments[0].range_start = start + 1
                feat.segments[0].range_end = end
            return True
        return False

    def find_feature(self, name: str) -> Optional[Feature]:
        """Find a feature by name."""
        return self._data.find_feature(name)

    def list_features(self) -> List[Dict]:
        """Return a list of feature summary dicts."""
        return [
            {
                "name": f.name,
                "type": f.type,
                "start": f.start,
                "end": f.end,
                "strand": f.strand,
                "color": f.color,
                "length": f.length,
            }
            for f in self._data.features
        ]

    # --- Save / Reload ---

    def save(self, path: str) -> str:
        """Save the SnapGene file to an explicit path.

        Args:
            path: Output file path. Must be specified explicitly.

        Raises:
            ValueError: If no path is provided.
        """
        if not path:
            raise ValueError(
                "Output path must be specified explicitly. "
                "Use editor.save('output.dna') instead of editor.save()."
            )
        write_snapgene(self._data, path)

    def save_overwrite(self) -> str:
        """Save back to the original source file (overwrites source).

        Use this only when you intentionally want to modify the source file.
        """
        if not self._filepath:
            raise ValueError(
                "No source file to overwrite. "
                "Use editor.save('output.dna') to specify a path."
            )
        write_snapgene(self._data, self._filepath)

    def reload(self) -> str:
        """Reload from the original file (discards unsaved changes)."""
        if self._filepath:
            self._data = read_snapgene(self._filepath)


# =============================================================================
# 批量 sgRNA 特征编辑（唯一公开的批量接口）
# =============================================================================

def batch_from_lists_simple(
    save_names: List[str],
    new_names: List[str],
    new_seqs: List[str],
    fluor_list: List[str],
    output_dir: str = "",
    fluor_files: Dict[str, str] = None,
) -> str:
    """
    纯配置驱动的批量编辑，根据荧光标识自动选择对应模板文件。

    Parameters
    ----------
    save_names : List[str]
        保存文件名列表（不含 .dna 后缀，不含路径），一一对应。
    new_names : List[str]
        新的特征名称列表，一一对应。
    new_seqs : List[str]
        新序列列表，一一对应。
        写 "FLUOR" 表示使用对应荧光的原始 sgRNA 序列。
    fluor_list : List[str]
        荧光标识列表（A/B/C/G），与 new_seqs 一一对应。
    output_dir : str, optional
        输出目录，默认为当前目录。
    fluor_files : Dict[str, str], optional
        荧光到模板文件路径的映射。
        如果不传，则使用 FLUORESCENT_FILES（向后兼容）。

    Example
    -------
    >>> batch_from_lists_simple(
    ...     save_names=["ZFP001", "ZFP002", "ZFP003"],
    ...     new_names=["sgKlf2-1", "sgKlf2-2", "sgKlf2-3"],
    ...     new_seqs=["FLUOR", "AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT"],
    ...     fluor_list=["A", "A", "A"],
    ...     output_dir=r"D:\\plasmid\\temp",
    ... )
    """
    if fluor_files is None:
        fluor_files = FLUORESCENT_FILES

    if not (
        len(save_names)
        == len(new_names)
        == len(new_seqs)
        == len(fluor_list)
    ):
        raise ValueError("四个列表的长度必须一致")

    if not output_dir:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)

    # 收集处理数据和原始序列（用于生成日志和引物）
    processed = []  # (save_name, fluor_upper, old_seq, new_seq_upper)

    for save_name, new_name, new_seq, fluor in zip(
        save_names, new_names, new_seqs, fluor_list
    ):
        fluor_upper = fluor.upper()
        src_path = fluor_files.get(fluor_upper, "")
        if not src_path or not os.path.exists(src_path):
            raise FileNotFoundError(
                f"荧光 '{fluor_upper}' 的模板文件不存在: {src_path}"
            )

        editor = SnapGeneEditor(src_path)

        # 移除无关特征
        for feat_name in ["source"]:
            f = editor.find_feature(feat_name)
            if f:
                editor.remove_feature(feat_name)

        feat = editor.find_feature("sgRNA")
        if not feat:
            raise ValueError(
                f"荧光 '{fluor_upper}' 的模板文件中未找到 'sgRNA' 特征"
            )

        # 解析序列
        if new_seq.strip().upper() == "FLUOR":
            new_seq_upper = editor.sequence[feat.start:feat.end]
        else:
            new_seq_upper = new_seq.upper()

        old_seq = editor.sequence[feat.start:feat.end]

        editor.replace_in_feature("sgRNA", new_seq_upper)
        editor.rename_feature("sgRNA", new_name)

        out_file = os.path.join(output_dir, f"{save_name}.dna")
        editor.save(out_file)

        processed.append((save_name, fluor_upper, old_seq, new_seq_upper))

    # 生成输出
    output_parts = []
    for save_name, fluor_upper, old_seq, new_seq_upper in processed:
        output_parts.append(
            f"[{save_name}] {fluor_upper}: {old_seq} -> {new_seq_upper}"
        )

    output_parts.append("要做这些克隆，请你复制这些序列粘贴到 primer list 中")
    output_parts.append("-" * 40)

    for save_name, fluor_upper, old_seq, new_seq_upper in processed:
        rc_seq = revcomp(new_seq_upper)
        output_parts.append(new_seq_upper)
        output_parts.append(rc_seq)

    result = "\n".join(output_parts)
    print(result)  # 在 Jupyter 中打印会正确换行
    return result
