"""Data models for SnapGene file contents."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class FeatureSegment:
    """A segment within a feature (for multi-segment features)."""

    range_start: int
    range_end: int
    color: str = "#999999"
    segment_type: str = "standard"
    name: str = ""


@dataclass
class Feature:
    """A feature/annotation on a DNA sequence."""

    name: str
    type: str = "misc_feature"
    start: int = 0
    end: int = 0
    strand: int = 0  # 1: forward, -1: reverse, 0: none
    color: str = "#999999"
    segments: List[FeatureSegment] = field(default_factory=list)
    qualifiers: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if not self.segments:
            self.segments = [
                FeatureSegment(range_start=self.start + 1, range_end=self.end, color=self.color)
            ]

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Primer:
    """A primer annotation."""

    name: str
    sequence: str
    start: int = 0
    end: int = 0
    type: str = "Primer"
    color: str = "#999999"


@dataclass
class Note:
    """A key-value note/annotation on the sequence."""

    key: str
    value: Any


@dataclass
class SnapGeneFile:
    """Represents the full contents of a SnapGene .dna file."""

    name: str = "Untitled"
    sequence: str = ""
    features: List[Feature] = field(default_factory=list)
    primers: List[Primer] = field(default_factory=list)
    notes: Dict[str, Any] = field(default_factory=dict)

    is_dna: bool = True
    export_version: int = 2
    import_version: int = 2

    # Sequence properties
    is_circular: bool = False
    is_double_stranded: bool = True
    dam_methylated: bool = False
    dcm_methylated: bool = False
    eco_ki_methylated: bool = False

    # Raw blocks we don't understand (preserved for round-trip)
    _raw_blocks: List[Dict[str, Any]] = field(default_factory=list)

    def add_feature(self, feature: Feature) -> None:
        """Add a feature to the sequence."""
        self.features.append(feature)

    def remove_feature(self, name: str) -> None:
        """Remove a feature by name."""
        self.features = [f for f in self.features if f.name != name]

    def find_feature(self, name: str) -> Optional[Feature]:
        """Find a feature by name."""
        for f in self.features:
            if f.name == name:
                return f
        return None

    @property
    def length(self) -> int:
        return len(self.sequence)

    def reverse_complement(self) -> "SnapGeneFile":
        """Return a new SnapGeneFile with the reverse complement sequence."""
        from snapgene_edit.utils import reverse_complement as revcomp

        new_seq = revcomp(self.sequence)
        new_features = []
        seq_len = len(self.sequence)
        for f in self.features:
            new_f = Feature(
                name=f.name,
                type=f.type,
                start=seq_len - f.end,
                end=seq_len - f.start,
                strand=-f.strand if f.strand != 0 else 0,
                color=f.color,
                qualifiers=f.qualifiers,
            )
            new_features.append(new_f)

        return SnapGeneFile(
            name=self.name,
            sequence=new_seq,
            features=new_features,
            primers=self.primers,
            notes=self.notes,
            is_circular=self.is_circular,
            is_double_stranded=self.is_double_stranded,
        )
