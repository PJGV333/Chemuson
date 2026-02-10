"""Regresiones para selección rectangular ceñida."""

import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import BondStereo, BondStyle
from gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_rect_selection_uses_item_centers_to_avoid_neighbor_leak():
    """Rectángulo de selección no debe arrastrar enlaces vecinos por solapamiento parcial."""
    canvas = ChemusonCanvas()

    a1 = canvas.model.add_atom("N", 120.0, 120.0).id
    a2 = canvas.model.add_atom("N", 160.0, 120.0).id
    a3 = canvas.model.add_atom("N", 174.0, 120.0).id
    a4 = canvas.model.add_atom("N", 214.0, 120.0).id
    b1 = canvas.model.add_bond(
        a1,
        a2,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    ).id
    b2 = canvas.model.add_bond(
        a3,
        a4,
        order=1,
        style=BondStyle.PLAIN,
        stereo=BondStereo.NONE,
        is_aromatic=False,
    ).id
    canvas._rebuild_items_from_model()

    # El rectángulo toca una porción del segundo enlace, pero su centro queda fuera.
    canvas._begin_selection_drag(QPointF(108.0, 108.0), free_select=False, additive=False)
    canvas._last_scene_pos = QPointF(171.0, 132.0)
    canvas._finalize_selection_drag()

    assert b1 in canvas.state.selected_bonds
    assert b2 not in canvas.state.selected_bonds
