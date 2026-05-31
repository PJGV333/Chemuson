"""Servicios de limpieza 2D desacoplados de la GUI."""

from .v2 import Clean2DParameters, optimize_clean2d_positions
from .safety import (
    Clean2DQualityReport,
    evaluate_clean2d_layout,
    is_clean2d_candidate_safe,
    has_cycles,
    count_new_bond_crossings,
    min_nonbonded_distance,
    ring_degeneracy_score,
    max_atom_displacement,
    bond_length_stats,
)
from .length_only import (
    length_only_polish,
    structure_preserving_geometry_polish,
    structure_preserving_length_polish,
)

__all__ = [
    "Clean2DParameters",
    "optimize_clean2d_positions",
    "Clean2DQualityReport",
    "evaluate_clean2d_layout",
    "is_clean2d_candidate_safe",
    "has_cycles",
    "count_new_bond_crossings",
    "min_nonbonded_distance",
    "ring_degeneracy_score",
    "max_atom_displacement",
    "bond_length_stats",
    "length_only_polish",
    "structure_preserving_geometry_polish",
    "structure_preserving_length_polish",
]
