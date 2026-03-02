"""Pruebas básicas de documentos por pestañas en la ventana principal."""

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication, QInputDialog

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio.persistence import PersistenceManager
from chemuson.chemio import rdkit_io
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.main_window import ChemusonWindow


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def test_new_creates_independent_tab() -> None:
    window = ChemusonWindow()
    initial = window.tabs.count()
    first_canvas = window.canvas

    window._on_file_new()

    assert window.tabs.tabsClosable()
    assert window.tabs.count() == initial + 1
    assert window.canvas is not first_canvas
    assert window.tabs.indexOf(first_canvas) >= 0


def test_tab_state_is_independent() -> None:
    window = ChemusonWindow()
    first = window.canvas
    first.state.show_implicit_carbons = True
    first.state.use_aromatic_circles = True
    first.set_numbering_mode("structures")
    first.set_numbering_enabled(True)

    window._on_file_new()
    second = window.canvas
    second.state.show_implicit_carbons = False
    second.state.use_aromatic_circles = False
    second.set_numbering_mode("atoms")
    second.set_numbering_enabled(False)

    window.tabs.setCurrentWidget(first)
    assert window.canvas is first
    assert window.canvas.state.show_implicit_carbons is True
    assert window.action_show_carbons.isChecked() is True
    assert window.action_aromatic_circles.isChecked() is True
    assert window.action_numbering_enabled.isChecked() is True

    window.tabs.setCurrentWidget(second)
    assert window.canvas is second
    assert window.canvas.state.show_implicit_carbons is False
    assert window.action_show_carbons.isChecked() is False
    assert window.action_aromatic_circles.isChecked() is False
    assert window.action_numbering_enabled.isChecked() is False


def test_open_file_path_creates_new_tab(tmp_path) -> None:
    sample_path = tmp_path / "sample.cmsn"
    source_canvas = ChemusonCanvas()
    PersistenceManager.save_to_file(str(sample_path), source_canvas)

    window = ChemusonWindow()
    initial = window.tabs.count()

    window._open_file_path(str(sample_path))

    assert window.tabs.count() == initial + 1
    assert window._canvas_file_paths[window.canvas] == os.path.abspath(str(sample_path))


def test_tab_change_clears_active_tool_selection() -> None:
    window = ChemusonWindow()
    first = window.canvas

    window._on_file_new()
    second = window.canvas

    window.tabs.setCurrentWidget(first)
    window.toolbar.bond_action.trigger()
    assert window.toolbar.action_group.checkedAction() is not None
    assert window.canvas.current_tool == "tool_bond"

    window.tabs.setCurrentWidget(second)

    assert window.canvas is second
    assert not any(action.isChecked() for action in window.toolbar.action_group.actions())
    assert window.canvas.current_tool == "tool_none"
    assert window.canvas.state.active_tool == "tool_none"


def test_import_smiles_inserts_without_clearing_and_undo(monkeypatch) -> None:
    window = ChemusonWindow()

    # Estructura previa en el documento (fuera de undo stack).
    a1 = window.canvas.model.add_atom("C", 10.0, 10.0)
    a2 = window.canvas.model.add_atom("C", 40.0, 10.0)
    b = window.canvas.model.add_bond(a1.id, a2.id, order=1)
    window.canvas.add_atom_item(a1)
    window.canvas.add_atom_item(a2)
    window.canvas.add_bond_item(b)

    baseline_atoms = set(window.canvas.model.atoms.keys())
    baseline_count = len(baseline_atoms)

    def _fake_get_text(*_args, **_kwargs):
        return "CC", True

    def _fake_smiles_to_molgraph(_smiles: str) -> MolGraph:
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)
        return graph

    monkeypatch.setattr(QInputDialog, "getText", staticmethod(_fake_get_text))
    monkeypatch.setattr(rdkit_io, "smiles_to_molgraph", _fake_smiles_to_molgraph)

    window._on_import_smiles()
    current_atoms = set(window.canvas.model.atoms.keys())
    assert len(current_atoms) == baseline_count + 2
    assert baseline_atoms.issubset(current_atoms)

    window._on_undo()
    assert set(window.canvas.model.atoms.keys()) == baseline_atoms


def test_insert_molgraph_skips_duplicate_bonds() -> None:
    canvas = ChemusonCanvas()
    source = MolGraph()
    a1 = source.add_atom("C", 0.0, 0.0)
    a2 = source.add_atom("C", 1.0, 0.0)
    source.add_bond(a1.id, a2.id, order=1)
    source.add_bond(a1.id, a2.id, order=1)

    canvas._insert_molgraph(source)
    assert len(canvas.model.bonds) == 1


def test_gallery_template_selection_enters_insert_mode() -> None:
    window = ChemusonWindow()
    try:
        baseline_atoms = len(window.canvas.model.atoms)

        templates = window.template_library.list_templates()
        benzene = next((tpl for tpl in templates if tpl.get("name") == "Benceno"), None)
        assert benzene is not None

        window._on_template_selected_from_gallery({"id": benzene["id"]})

        assert len(window.canvas.model.atoms) == baseline_atoms
        assert window.canvas._pending_template_graph is not None
        assert window.canvas._pending_template_label
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()
