"""Pruebas del panel interactivo de validación química."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.main_window import ChemusonWindow


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
        assert window.validation_dock.issue_table.item(0, 0).text() == str(atom_id)
        assert "VALENCE" in window.validation_dock.issue_table.item(0, 2).text()
        assert window.validation_dock.issue_table.item(0, 4).text()
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
