"""Regresiones para tamaño de etiquetas según selección estructural."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_apply_label_font_size_uses_selected_atoms_only():
    """El tamaño exacto debe aplicarse solo a etiquetas de átomos seleccionados."""
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("N", 100.0, 100.0)
    a2 = canvas.model.add_atom("O", 140.0, 100.0)
    a3 = canvas.model.add_atom("Cl", 180.0, 100.0)
    canvas._rebuild_items_from_model()

    canvas.atom_items[a1.id].setSelected(True)
    canvas.atom_items[a2.id].setSelected(True)
    canvas._sync_selection_from_scene()

    assert canvas.apply_label_font_size(20.0)

    assert canvas.state.label_font_size == pytest.approx(11.0)
    assert canvas._effective_label_font_size(a1.id) == pytest.approx(20.0)
    assert canvas._effective_label_font_size(a2.id) == pytest.approx(20.0)
    assert canvas._effective_label_font_size(a3.id) == pytest.approx(11.0)
    assert canvas.model.get_atom(a1.id).label_scale == pytest.approx(20.0 / 11.0)
    assert canvas.model.get_atom(a2.id).label_scale == pytest.approx(20.0 / 11.0)
    assert canvas.model.get_atom(a3.id).label_scale is None

    canvas.undo_stack.undo()

    assert canvas._effective_label_font_size(a1.id) == pytest.approx(11.0)
    assert canvas._effective_label_font_size(a2.id) == pytest.approx(11.0)
    assert canvas.model.get_atom(a1.id).label_scale is None
    assert canvas.model.get_atom(a2.id).label_scale is None


def test_adjust_label_font_size_uses_selected_bond_endpoints_only():
    """Aumentar/reducir tamaño debe afectar solo la estructura seleccionada."""
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("N", 100.0, 100.0)
    a2 = canvas.model.add_atom("O", 140.0, 100.0)
    a3 = canvas.model.add_atom("Cl", 180.0, 100.0)
    bond = canvas.model.add_bond(a1.id, a2.id)
    canvas._rebuild_items_from_model()

    canvas.bond_items[bond.id].setSelected(True)
    canvas._sync_selection_from_scene()

    assert canvas.adjust_label_font_size(2.0)

    assert canvas.state.label_font_size == pytest.approx(11.0)
    assert canvas._effective_label_font_size(a1.id) == pytest.approx(13.0)
    assert canvas._effective_label_font_size(a2.id) == pytest.approx(13.0)
    assert canvas._effective_label_font_size(a3.id) == pytest.approx(11.0)


def test_apply_label_font_size_without_selection_remains_global():
    """Sin selección estructural, el cambio de tamaño debe seguir siendo global."""
    canvas = ChemusonCanvas()
    a1 = canvas.model.add_atom("N", 100.0, 100.0)
    a2 = canvas.model.add_atom("O", 140.0, 100.0)
    canvas._rebuild_items_from_model()

    assert canvas.apply_label_font_size(18.0)

    assert canvas.state.label_font_size == pytest.approx(18.0)
    assert canvas._effective_label_font_size(a1.id) == pytest.approx(18.0)
    assert canvas._effective_label_font_size(a2.id) == pytest.approx(18.0)
