"""Regresiones de GUI para el menú de rotación de fragmentos con pivote."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.main_window import ChemusonWindow


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def test_fragment_pivot_menu_defines_rotates_and_clears() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom_a = canvas.model.add_atom("C", 80.0, 100.0)
        atom_b = canvas.model.add_atom("N", 120.0, 100.0)
        atom_c = canvas.model.add_atom("C", 160.0, 100.0)
        atom_d = canvas.model.add_atom("C", 200.0, 100.0)

        canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
        canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
        canvas.model.add_bond(atom_c.id, atom_d.id, order=1)
        canvas._rebuild_items_from_model()

        canvas.atom_items[atom_b.id].setSelected(True)
        canvas._sync_selection_from_scene()

        assert window.action_fragment_pivot_set.isEnabled() is True
        assert f"átomo {atom_b.id}" in window.action_fragment_pivot_set.text()

        window._on_set_fragment_pivot()

        assert canvas.fragment_pivot_atom_id() == atom_b.id
        assert f"átomo {atom_b.id}" in window.fragment_rotate_menu.title()
        assert f"átomo {atom_b.id}" in window.action_fragment_pivot_clear.text()

        canvas.scene.clearSelection()
        canvas.atom_items[atom_c.id].setSelected(True)
        canvas.atom_items[atom_d.id].setSelected(True)
        canvas._sync_selection_from_scene()

        assert window.action_fragment_rotate_plus.isEnabled() is True

        window._on_rotate_fragment(90.0)

        assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
            (120.0, 140.0)
        )
        assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
            (120.0, 180.0)
        )

        window._on_undo()

        assert (canvas.model.get_atom(atom_c.id).x, canvas.model.get_atom(atom_c.id).y) == pytest.approx(
            (160.0, 100.0)
        )
        assert (canvas.model.get_atom(atom_d.id).x, canvas.model.get_atom(atom_d.id).y) == pytest.approx(
            (200.0, 100.0)
        )

        window._on_clear_fragment_pivot()

        assert canvas.fragment_pivot_atom_id() is None
        assert window.action_fragment_pivot_clear.isEnabled() is False
        assert window.fragment_rotate_menu.title() == "Fragmento con pivote"
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()


def test_fragment_pivot_menu_enablement_tracks_selection() -> None:
    window = ChemusonWindow()
    try:
        canvas = window.canvas
        atom_a = canvas.model.add_atom("C", 80.0, 100.0)
        atom_b = canvas.model.add_atom("C", 120.0, 100.0)
        atom_c = canvas.model.add_atom("C", 160.0, 100.0)
        bond_ab = canvas.model.add_bond(atom_a.id, atom_b.id, order=1)
        canvas.model.add_bond(atom_b.id, atom_c.id, order=1)
        canvas._rebuild_items_from_model()

        assert window.action_fragment_pivot_set.isEnabled() is False
        assert window.action_fragment_rotate_plus.isEnabled() is False

        canvas.atom_items[atom_a.id].setSelected(True)
        canvas.atom_items[atom_b.id].setSelected(True)
        canvas._sync_selection_from_scene()

        assert window.action_fragment_pivot_set.isEnabled() is False

        canvas.scene.clearSelection()
        canvas.atom_items[atom_b.id].setSelected(True)
        canvas._sync_selection_from_scene()

        assert window.action_fragment_pivot_set.isEnabled() is True

        window._on_set_fragment_pivot()

        canvas.scene.clearSelection()
        canvas._sync_selection_from_scene()

        assert window.action_fragment_rotate_plus.isEnabled() is False

        canvas.bond_items[bond_ab.id].setSelected(True)
        canvas._sync_selection_from_scene()

        assert window.action_fragment_rotate_plus.isEnabled() is True
        assert window.action_fragment_rotate_minus.isEnabled() is True
        assert window.action_fragment_invert.isEnabled() is True
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()
