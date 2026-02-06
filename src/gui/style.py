"""
Drawing style presets for Chemuson.
"""
from __future__ import annotations

from dataclasses import dataclass
from PyQt6.QtCore import Qt


@dataclass(frozen=True)
class DrawingStyle:
    bond_length_px: float
    stroke_px: float
    bond_color: str
    atom_fill_color: str
    atom_stroke_color: str
    double_offset_px: float
    inner_trim_px: float
    bond_end_shrink_px: float
    wedge_width_px: float
    hash_count: int
    hash_min_px: float
    hash_max_px: float
    hash_stroke_px: float
    cap_style: Qt.PenCapStyle
    join_style: Qt.PenJoinStyle


CHEMDOODLE_LIKE = DrawingStyle(
    bond_length_px=40.0,
    stroke_px=2.6,
    bond_color="#000000",
    atom_fill_color="#FFFFFF",
    atom_stroke_color="#333333",
    double_offset_px=4.2,
    inner_trim_px=7.0,
    bond_end_shrink_px=0.0,
    wedge_width_px=12.0,
    hash_count=6,
    hash_min_px=2.0,
    hash_max_px=10.0,
    hash_stroke_px=2.2,
    cap_style=Qt.PenCapStyle.RoundCap,
    join_style=Qt.PenJoinStyle.RoundJoin,
)
