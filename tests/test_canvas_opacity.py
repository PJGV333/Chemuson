"""Pruebas para transparencia/opacidad de selección y canvas completo."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import TextAnnotationItem
from chemuson.gui.main_window import ChemusonWindow


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def test_apply_opacity_uses_selection_then_canvas_default():
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("N", 100.0, 100.0)
    atom_b = canvas.model.add_atom("O", 140.0, 100.0)
    bond = canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    text_item = TextAnnotationItem("Nota", 200.0, 100.0)
    canvas.add_text_item(text_item)
    arrow = canvas.add_arrow_item(QPointF(220.0, 100.0), QPointF(280.0, 100.0), "forward")

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas._sync_selection_from_scene()

    assert canvas.apply_opacity_percent(35)
    assert canvas.model.get_atom(atom_a.id).opacity == pytest.approx(0.35)
    assert canvas.effective_atom_opacity(atom_a.id) == pytest.approx(0.35)
    assert canvas.effective_atom_opacity(atom_b.id) == pytest.approx(1.0)
    assert canvas.effective_bond_opacity(bond.id) == pytest.approx(1.0)
    assert canvas.effective_item_opacity(text_item) == pytest.approx(1.0)
    assert canvas.effective_item_opacity(arrow) == pytest.approx(1.0)

    canvas.scene.clearSelection()
    canvas._sync_selection_from_scene()

    assert canvas.current_opacity_percent() == 100
    assert canvas.apply_opacity_percent(50)
    assert canvas.canvas_default_opacity() == pytest.approx(0.5)
    assert canvas.model.get_atom(atom_a.id).opacity is None
    assert canvas.model.get_atom(atom_b.id).opacity is None
    assert canvas.model.get_bond(bond.id).opacity is None
    assert canvas.effective_atom_opacity(atom_a.id) == pytest.approx(0.5)
    assert canvas.effective_atom_opacity(atom_b.id) == pytest.approx(0.5)
    assert canvas.effective_bond_opacity(bond.id) == pytest.approx(0.5)
    assert canvas.effective_item_opacity(text_item) == pytest.approx(0.5)
    assert canvas.effective_item_opacity(arrow) == pytest.approx(0.5)

    inherited_arrow = canvas.add_arrow_item(QPointF(300.0, 100.0), QPointF(360.0, 100.0), "forward")
    assert canvas.effective_item_opacity(inherited_arrow) == pytest.approx(0.5)

    canvas.undo_stack.undo()

    assert canvas.canvas_default_opacity() == pytest.approx(1.0)
    assert canvas.model.get_atom(atom_a.id).opacity == pytest.approx(0.35)
    assert canvas.effective_atom_opacity(atom_a.id) == pytest.approx(0.35)
    assert canvas.effective_atom_opacity(atom_b.id) == pytest.approx(1.0)
    assert canvas.effective_bond_opacity(bond.id) == pytest.approx(1.0)
    assert canvas.effective_item_opacity(text_item) == pytest.approx(1.0)


def test_opacity_persists_and_selection_payload_preserves_visual_result():
    canvas = ChemusonCanvas()
    assert canvas.apply_opacity_percent(40)

    atom_a = canvas.model.add_atom("C", 60.0, 80.0)
    atom_b = canvas.model.add_atom("C", 100.0, 80.0)
    bond = canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    text_item = TextAnnotationItem("Bloque", 150.0, 70.0)
    canvas.add_text_item(text_item)
    arrow = canvas.add_arrow_item(QPointF(180.0, 80.0), QPointF(240.0, 80.0), "forward")

    canvas.bond_items[bond.id].setSelected(True)
    text_item.setSelected(True)
    arrow.setSelected(True)
    canvas._sync_selection_from_scene()

    assert canvas.apply_opacity_percent(80)

    payload = canvas._build_selection_payload()
    assert payload is not None
    assert payload["bonds"][0]["opacity"] == pytest.approx(0.8)
    assert payload["texts"][0]["opacity"] == pytest.approx(0.8)
    assert payload["arrows"][0]["opacity"] == pytest.approx(0.8)

    pasted = ChemusonCanvas()
    pasted._last_scene_pos = QPointF(300.0, 220.0)
    pasted._paste_selection_payload(payload)

    pasted_bond = next(iter(pasted.model.bonds.values()))
    pasted_text = pasted._all_scalable_text_items()[0]
    assert pasted.effective_bond_opacity(pasted_bond.id) == pytest.approx(0.8)
    assert pasted.effective_item_opacity(pasted.arrow_items[0]) == pytest.approx(0.8)
    assert pasted.effective_item_opacity(pasted_text) == pytest.approx(0.8)

    data = PersistenceManager.save_to_dict(canvas)
    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(data, restored)

    restored_atom_a = restored.model.get_atom(atom_a.id)
    restored_atom_b = restored.model.get_atom(atom_b.id)
    restored_bond = next(iter(restored.model.bonds.values()))
    restored_text = restored._all_scalable_text_items()[0]

    assert restored.canvas_default_opacity() == pytest.approx(0.4)
    assert restored_atom_a.opacity is None
    assert restored_atom_b.opacity is None
    assert restored.effective_atom_opacity(restored_atom_a.id) == pytest.approx(0.4)
    assert restored.effective_atom_opacity(restored_atom_b.id) == pytest.approx(0.4)
    assert restored_bond.opacity == pytest.approx(0.8)
    assert restored.effective_bond_opacity(restored_bond.id) == pytest.approx(0.8)
    assert restored.effective_item_opacity(restored.arrow_items[0]) == pytest.approx(0.8)
    assert restored.effective_item_opacity(restored_text) == pytest.approx(0.8)


def test_toolbar_opacity_control_routes_to_selection_or_canvas():
    window = ChemusonWindow()
    try:
        atom_a = window.canvas.model.add_atom("N", 80.0, 80.0)
        atom_b = window.canvas.model.add_atom("N", 120.0, 80.0)
        window.canvas.model.add_bond(atom_a.id, atom_b.id)
        window.canvas._rebuild_items_from_model()

        window.canvas.atom_items[atom_a.id].setSelected(True)
        window.canvas._sync_selection_from_scene()
        window.text_toolbar.opacity_spin.setValue(30)

        assert window.canvas.effective_atom_opacity(atom_a.id) == pytest.approx(0.3)
        assert window.canvas.effective_atom_opacity(atom_b.id) == pytest.approx(1.0)

        window.canvas.scene.clearSelection()
        window.canvas._sync_selection_from_scene()
        window.text_toolbar.opacity_spin.setValue(55)

        assert window.canvas.canvas_default_opacity() == pytest.approx(0.55)
        assert window.canvas.effective_atom_opacity(atom_b.id) == pytest.approx(0.55)
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()
