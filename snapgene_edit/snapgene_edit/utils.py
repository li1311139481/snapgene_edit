"""Utility functions for DNA sequence manipulation."""

_COMPLEMENT_TABLE = {
    "A": "T", "B": "V", "C": "G", "D": "H",
    "G": "C", "H": "D", "K": "M", "M": "K",
    "N": "N", "R": "Y", "S": "S", "T": "A",
    "U": "A", "V": "B", "W": "W", "X": "X",
    "Y": "R",
    "a": "t", "b": "v", "c": "g", "d": "h",
    "g": "c", "h": "d", "k": "m", "m": "k",
    "n": "n", "r": "y", "s": "s", "t": "a",
    "u": "a", "v": "b", "w": "w", "x": "x",
    "y": "r",
}


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence."""
    return "".join(_COMPLEMENT_TABLE.get(c, c) for c in seq)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return complement(seq)[::-1]


def guess_type(seq: str) -> str:
    """Infer the type of a sequence."""
    import re
    if re.match(r"^[atgcn.]+$", seq, re.IGNORECASE):
        return "dna"
    elif re.match(r"^[augcn.]+$", seq, re.IGNORECASE):
        return "rna"
    return "unknown"


def hex_color_to_int(color: str) -> int:
    """Convert hex color string like '#FF0000' to integer."""
    color = color.lstrip("#")
    if len(color) == 6:
        r, g, b = color[0:2], color[2:4], color[4:6]
        return (int(r, 16) << 16) | (int(g, 16) << 8) | int(b, 16)
    return 0


def int_to_hex_color(value: int) -> str:
    """Convert integer to hex color string."""
    return f"#{value:06X}"
