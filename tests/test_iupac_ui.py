"""Pruebas UI para el campo de nombre IUPAC y preferencias relacionadas."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication, QLabel

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import ChemState, MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.dialogs import PreferencesDialog
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.style import CHEMDOODLE_LIKE


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _ethane_graph() -> MolGraph:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 1.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    return graph


def test_canvas_analysis_uses_nombre_iupac_label() -> None:
    canvas = ChemusonCanvas()
    text = canvas._analysis_build_text(_ethane_graph(), "name")
    assert "Nombre IUPAC:" in text


def test_preferences_dialog_emits_naming_options() -> None:
    dialog = PreferencesDialog(
        ChemState(),
        CHEMDOODLE_LIKE,
        update_settings={},
        naming_settings={"advanced_enabled": True, "rdkit_isolated": True},
    )
    emitted = {}

    def _capture(payload: dict) -> None:
        emitted.update(payload)

    dialog.preferences_changed.connect(_capture)
    dialog.advanced_name_checkbox.setChecked(False)
    dialog.rdkit_isolated_checkbox.setChecked(False)
    dialog._on_accept()

    assert emitted.get("name_advanced_enabled") is False
    assert emitted.get("name_rdkit_isolated") is False


def test_preferences_dialog_explains_update_delivery_modes() -> None:
    dialog = PreferencesDialog(
        ChemState(),
        CHEMDOODLE_LIKE,
        update_settings={},
        naming_settings={"advanced_enabled": True, "rdkit_isolated": True},
    )

    all_text = "\n".join(label.text() for label in dialog.findChildren(QLabel))
    assert "Buscar actualizaciones cuando esta edición lo permita" == dialog.update_enabled_checkbox.text()
    assert "Notificar antes de aplicar" == dialog.update_mode_combo.itemText(0)
    assert "Preparar en silencio cuando sea posible" == dialog.update_mode_combo.itemText(1)
    assert "auto-update real" in all_text.lower()
    assert "flatpak update io.github.pjgv333.chemuson" in all_text.lower()


def test_main_window_shows_iupac_status_field() -> None:
    window = ChemusonWindow()
    assert window._iupac_name_label.text().startswith("Nombre IUPAC:")
