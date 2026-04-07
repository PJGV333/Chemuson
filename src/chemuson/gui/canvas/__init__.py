"""API publica estable del paquete canvas."""
from __future__ import annotations

from PyQt6.QtWidgets import QInputDialog

from chemuson.chemio.rdkit_io import molgraph_to_molfile, molgraph_to_smiles

from .canvas_constants import (
    AROMATIC_CIRCLE_ATOMS_ROLE,
    BRANCH_ROTATION_NOOP_TOLERANCE_DEG,
    BRANCH_ROTATION_STEP_DEG,
    FRAGMENT_ROTATION_STEP_DEG,
)
from .canvas_view import ChemusonCanvas

__all__ = [
    "AROMATIC_CIRCLE_ATOMS_ROLE",
    "BRANCH_ROTATION_NOOP_TOLERANCE_DEG",
    "BRANCH_ROTATION_STEP_DEG",
    "ChemusonCanvas",
    "FRAGMENT_ROTATION_STEP_DEG",
    "molgraph_to_molfile",
    "molgraph_to_smiles",
    "QInputDialog",
]
