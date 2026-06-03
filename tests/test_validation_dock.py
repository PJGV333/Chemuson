"""Pruebas del panel interactivo de validación química."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.controllers import ValidationController


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _add_overvalent_carbon(window: ChemusonWindow) -> int:
    carbon = window.canvas.model.add_atom("C", 0.0, 0.0)
    for idx in range(5):
        hydrogen = window.canvas.model.add_atom("H", float(idx + 1) * 40.0, 0.0)
        window.canvas.model.add_bond(carbon.id, hydrogen.id, order=1)
    window.canvas._rebuild_items_from_model()
    return carbon.id


def test_validate_structure_populates_validation_dock() -> None:
    window = ChemusonWindow()
    try:
        atom_id = _add_overvalent_carbon(window)

        window._on_validate_structure()

        assert not window.validation_dock.isHidden()
        assert window.validation_dock.issue_table.rowCount() == 1
        assert window.validation_dock.issue_table.item(0, 0).text() == f"Átomo {atom_id} (C)"
        assert "VALENCE" in window.validation_dock.issue_table.item(0, 2).text()
        assert "enlaces" in window.validation_dock.issue_table.item(0, 3).text()
        assert window.validation_dock.issue_table.item(0, 5).text()
    finally:
        window.close()


def test_validation_dock_selection_selects_canvas_atom() -> None:
    window = ChemusonWindow()
    try:
        atom_id = _add_overvalent_carbon(window)
        window._on_validate_structure()

        window.validation_dock.issue_table.selectRow(0)
        QApplication.processEvents()

        assert window.canvas.atom_items[atom_id].isSelected()
    finally:
        window.close()


def test_validation_dock_copies_and_exports_report(tmp_path) -> None:
    window = ChemusonWindow()
    try:
        _add_overvalent_carbon(window)
        window._on_validate_structure()

        text = window.validation_dock.copy_report_tsv()
        assert "observed_valence" in text
        assert "VALENCE_EXCEEDED" in text
        assert QApplication.clipboard().text() == text

        out = tmp_path / "validation.json"
        window.validation_dock.export_report(str(out))
        assert "VALENCE_EXCEEDED" in out.read_text(encoding="utf-8")
    finally:
        window.close()


def test_validation_controller_adjusts_charge_undoably() -> None:
    window = ChemusonWindow()
    try:
        n = window.canvas.model.add_atom("N", 0.0, 0.0)
        for idx in range(4):
            c = window.canvas.model.add_atom("C", float(idx + 1) * 40.0, 0.0)
            window.canvas.model.add_bond(n.id, c.id, order=1)
        window.canvas._rebuild_items_from_model()

        window._on_validate_structure()
        window.validation_dock.issue_table.selectRow(0)
        QApplication.processEvents()
        labels = [window.validation_dock.action_combo.itemText(i) for i in range(window.validation_dock.action_combo.count())]
        assert any("carga formal" in label for label in labels)

        window.validation_dock.btn_apply_correction.click()
        QApplication.processEvents()

        assert window.canvas.model.get_atom(n.id).charge == 1
        assert n.id not in window.canvas.model.validate()
        assert window.canvas.undo_stack.count() == 1
        window.canvas.undo_stack.undo()
        assert window.canvas.model.get_atom(n.id).charge == 0
    finally:
        window.close()


def test_validation_controller_reduces_selected_bond_when_safe() -> None:
    window = ChemusonWindow()
    try:
        c = window.canvas.model.add_atom("C", 0.0, 0.0)
        o = window.canvas.model.add_atom("O", 40.0, 0.0)
        selected = window.canvas.model.add_bond(c.id, o.id, order=2)
        for idx in range(3):
            h = window.canvas.model.add_atom("H", 0.0, float(idx + 1) * 40.0)
            window.canvas.model.add_bond(c.id, h.id, order=1)
        window.canvas._rebuild_items_from_model()
        window.canvas.state.selected_bonds = {selected.id}
        window.canvas.validate_structure()
        issue = window.canvas.current_validation_issues()[c.id]

        actions = ValidationController.available_actions(window.canvas, issue)
        assert any(action["id"] == "reduce_selected_bond" for action in actions)
        result = ValidationController.apply_correction(window.canvas, c.id, "reduce_selected_bond")

        assert result.applied
        assert window.canvas.model.get_bond(selected.id).order == 1
        assert c.id not in window.canvas.model.validate()
        window.canvas.undo_stack.undo()
        assert window.canvas.model.get_bond(selected.id).order == 2
    finally:
        window.close()


def test_validation_controller_clears_assigned_h_when_safe() -> None:
    window = ChemusonWindow()
    try:
        o = window.canvas.model.add_atom("O", 0.0, 0.0, explicit_h=1, is_explicit=True)
        c1 = window.canvas.model.add_atom("C", 40.0, 0.0)
        c2 = window.canvas.model.add_atom("C", -40.0, 0.0)
        window.canvas.model.add_bond(o.id, c1.id, order=1)
        window.canvas.model.add_bond(o.id, c2.id, order=1)
        window.canvas._rebuild_items_from_model()
        window.canvas.validate_structure()
        issue = window.canvas.current_validation_issues()[o.id]

        actions = ValidationController.available_actions(window.canvas, issue)
        assert any(action["id"] == "clear_assigned_h" for action in actions)
        result = ValidationController.apply_correction(window.canvas, o.id, "clear_assigned_h")

        assert result.applied
        assert window.canvas.model.get_atom(o.id).explicit_h == 0
        assert o.id not in window.canvas.model.validate()
        window.canvas.undo_stack.undo()
        assert window.canvas.model.get_atom(o.id).explicit_h == 1
    finally:
        window.close()


def test_validation_navigation_syncs_dock_selection() -> None:
    window = ChemusonWindow()
    try:
        atom_id = _add_overvalent_carbon(window)

        window._on_navigate_validation_issue(1)

        selected = window.validation_dock.issue_table.selectedItems()
        assert selected
        assert int(selected[0].data(Qt.ItemDataRole.UserRole)) == atom_id
    finally:
        window.close()
