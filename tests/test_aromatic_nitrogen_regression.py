"""Regresiones para N aromático en el editor."""

import math

import pytest
from PyQt6.QtWidgets import QApplication


from chemuson.core.model import BondStereo, BondStyle
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import ChangeAtomCommand


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _build_aromatic_hexagon(
    canvas: ChemusonCanvas,
    center_x: float = 220.0,
    center_y: float = 220.0,
    radius: float = 42.0,
) -> list[int]:
    atom_ids: list[int] = []
    for idx in range(6):
        theta = math.radians(60.0 * idx)
        atom = canvas.model.add_atom(
            "C",
            center_x + radius * math.cos(theta),
            center_y + radius * math.sin(theta),
        )
        atom_ids.append(atom.id)
    for idx in range(6):
        canvas.model.add_bond(
            atom_ids[idx],
            atom_ids[(idx + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    return atom_ids


def test_change_aromatic_carbon_to_nitrogen_drops_spurious_h():
    """Sustituir un carbono aromático por N debe dibujar una piridina, no NH."""
    canvas = ChemusonCanvas()
    ring_atom_ids = _build_aromatic_hexagon(canvas)
    canvas._rebuild_items_from_model()

    target_atom_id = ring_atom_ids[0]
    canvas.undo_stack.push(ChangeAtomCommand(canvas.model, canvas, target_atom_id, "N"))

    assert canvas.model.implicit_h_count(target_atom_id) == 0
    assert canvas.atom_items[target_atom_id].label.toPlainText() == "N"
