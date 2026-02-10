"""Regresiones de visualización para círculos aromáticos."""

import math
import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStyle, BondStereo
from gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _build_hexagon(canvas: ChemusonCanvas, center_x: float = 220.0, center_y: float = 220.0):
    ids = []
    radius = 42.0
    for i in range(6):
        theta = math.radians(60.0 * i)
        atom = canvas.model.add_atom(
            "C",
            center_x + radius * math.cos(theta),
            center_y + radius * math.sin(theta),
        )
        ids.append(atom.id)
    return ids


def test_aromatic_circle_survives_single_wedge_styling():
    """Un anillo aromático con una cuña debe conservar el círculo aromático."""
    canvas = ChemusonCanvas()
    atom_ids = _build_hexagon(canvas)

    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()

    # Simula un caso real donde el enlace estilizado queda no-aromático.
    wedge_bond = canvas.model.get_bond(1)
    canvas.model.update_bond(
        wedge_bond.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 1


def test_aromatic_circle_survives_two_stereo_styled_ring_bonds():
    """El círculo aromático se mantiene con dos enlaces de anillo en cuña/hash."""
    canvas = ChemusonCanvas()
    atom_ids = _build_hexagon(canvas, center_x=340.0, center_y=220.0)

    for i in range(6):
        canvas.model.add_bond(
            atom_ids[i],
            atom_ids[(i + 1) % 6],
            order=1,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=True,
        )
    canvas._kekulize_aromatic_bonds()

    bond_1 = canvas.model.get_bond(1)
    bond_2 = canvas.model.get_bond(2)
    canvas.model.update_bond(
        bond_1.id,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        is_aromatic=False,
    )
    canvas.model.update_bond(
        bond_2.id,
        style=BondStyle.HASHED,
        stereo=BondStereo.DOWN,
        is_aromatic=False,
    )

    canvas._rebuild_items_from_model()
    canvas.state.use_aromatic_circles = True
    canvas.refresh_aromatic_circles()

    assert len(canvas.aromatic_circles) == 1
