"""Pruebas básicas de documentos por pestañas en la ventana principal."""

import os
import time

import pytest
from PyQt6.QtCore import QObject, Qt
from PyQt6.QtGui import QKeySequence
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QInputDialog, QMessageBox, QTabWidget, QTextEdit


from chemuson.core.model import MolGraph
from chemuson.chemio.persistence import PersistenceManager
from chemuson.chemio import rdkit_io
from chemuson.chemio import rdkit_safe
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.tab_manager import CanvasTabManager, _QtAutosaveManager
from chemuson.utils.autosave import AutosaveManager, AutosaveSerializer, AutosaveUndoStack


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


class FakeAutosaveController:
    def __init__(self) -> None:
        self.start_calls = 0
        self.stop_calls = 0
        self.original_paths: list[str | None] = []
        self.restart_calls = 0
        self.cancel_calls = 0
        self.cleanup_calls = 0

    def start(self) -> None:
        self.start_calls += 1

    def stop(self) -> None:
        self.stop_calls += 1

    def set_original_path(self, filepath: str | None) -> None:
        self.original_paths.append(filepath)

    def restart_debounce(self) -> None:
        self.restart_calls += 1

    def cancel_debounce(self) -> None:
        self.cancel_calls += 1

    def cleanup_after_save(self) -> None:
        self.cleanup_calls += 1


def test_tab_manager_injects_typed_autosave_collaborators() -> None:
    captured: list[tuple[QObject, ChemusonCanvas, AutosaveUndoStack, AutosaveSerializer[ChemusonCanvas]]] = []
    controller = FakeAutosaveController()

    def factory(
        parent: QObject,
        document: ChemusonCanvas,
        undo_stack: AutosaveUndoStack,
        serializer: AutosaveSerializer[ChemusonCanvas],
    ) -> FakeAutosaveController:
        captured.append((parent, document, undo_stack, serializer))
        return controller

    tabs = QTabWidget()
    parent = QObject()
    manager = CanvasTabManager(tabs, autosave_parent=parent, autosave_factory=factory)
    canvas = manager.create_document_tab()

    assert captured == [(parent, canvas, canvas.undo_stack, PersistenceManager.save_to_dict)]
    assert controller.start_calls == 1
    assert manager.autosave_managers[canvas] is controller
    assert manager.autosave_manager_for(canvas) is controller


def test_qt_autosave_adapter_owns_preserves_concrete_timer_configuration() -> None:
    window = ChemusonWindow()
    try:
        manager = window._tab_manager.autosave_manager_for(window.canvas)

        assert isinstance(manager, _QtAutosaveManager)
        assert isinstance(manager.core, AutosaveManager)
        assert manager.core._periodic_timer.parent() is manager
        assert manager.core._periodic_timer.interval() == AutosaveManager.AUTOSAVE_INTERVAL_MS
        assert manager.core._debounce_timer.isSingleShot()
        assert manager.core._debounce_timer.interval() == AutosaveManager.DEBOUNCE_INTERVAL_MS
        assert "__getattr__" not in _QtAutosaveManager.__dict__
    finally:
        window.close()


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


def test_copy_as_smiles_uses_selected_structure_only(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        first_a = window.canvas.model.add_atom("C", 10.0, 10.0)
        first_b = window.canvas.model.add_atom("C", 40.0, 10.0)
        second_a = window.canvas.model.add_atom("N", 110.0, 10.0)
        second_b = window.canvas.model.add_atom("N", 140.0, 10.0)
        window.canvas.model.add_bond(first_a.id, first_b.id, order=1)
        window.canvas.model.add_bond(second_a.id, second_b.id, order=1)
        window.canvas._rebuild_items_from_model()

        window.canvas.scene.clearSelection()
        window.canvas.atom_items[first_a.id].setSelected(True)
        window.canvas.atom_items[first_b.id].setSelected(True)
        window.canvas._sync_selection_from_scene()

        captured: dict[str, set[int]] = {}

        def _fake_molgraph_to_smiles(graph: MolGraph, timeout_s: float = 5.0) -> tuple[str, None]:
            captured["atom_ids"] = set(graph.atoms.keys())
            return "selected-smiles", None

        monkeypatch.setattr(rdkit_safe, "molgraph_to_smiles_isolated", _fake_molgraph_to_smiles)

        window._on_copy_as("smiles")

        assert captured["atom_ids"] == {first_a.id, first_b.id}
        assert second_a.id not in captured["atom_ids"]
        assert second_b.id not in captured["atom_ids"]
        assert QApplication.clipboard().text() == "selected-smiles"
    finally:
        window.close()


def test_export_smiles_uses_selected_structure_only(monkeypatch) -> None:
    window = ChemusonWindow()
    try:
        first_a = window.canvas.model.add_atom("C", 10.0, 10.0)
        first_b = window.canvas.model.add_atom("C", 40.0, 10.0)
        second_a = window.canvas.model.add_atom("N", 110.0, 10.0)
        second_b = window.canvas.model.add_atom("N", 140.0, 10.0)
        window.canvas.model.add_bond(first_a.id, first_b.id, order=1)
        window.canvas.model.add_bond(second_a.id, second_b.id, order=1)
        window.canvas._rebuild_items_from_model()

        window.canvas.scene.clearSelection()
        window.canvas.atom_items[first_a.id].setSelected(True)
        window.canvas.atom_items[first_b.id].setSelected(True)
        window.canvas._sync_selection_from_scene()

        captured: dict[str, object] = {}

        def _fake_molgraph_to_smiles(graph: MolGraph) -> str:
            captured["atom_ids"] = set(graph.atoms.keys())
            return "selected-smiles"

        def _fake_information(_parent, title: str, text: str) -> None:
            captured["dialog"] = (title, text)

        monkeypatch.setattr(
            rdkit_io,
            "molgraph_to_smiles_isolated_or_error",
            _fake_molgraph_to_smiles,
        )
        monkeypatch.setattr(QMessageBox, "information", _fake_information)

        window._on_export_smiles()
        deadline = time.monotonic() + 2.0
        while "dialog" not in captured and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)

        assert captured["atom_ids"] == {first_a.id, first_b.id}
        assert second_a.id not in captured["atom_ids"]
        assert second_b.id not in captured["atom_ids"]
        assert captured["dialog"] == ("SMILES", "selected-smiles")
    finally:
        window.close()


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


def test_external_rich_text_editor_uses_standard_text_toolbar() -> None:
    window = ChemusonWindow()
    try:
        editor = QTextEdit()
        editor.setPlainText("O2")
        window._set_external_text_editor(editor)

        cursor = editor.textCursor()
        cursor.setPosition(1)
        cursor.movePosition(cursor.MoveOperation.End, cursor.MoveMode.KeepAnchor)
        editor.setTextCursor(cursor)

        window._on_text_format_changed(
            "Arial",
            10,
            False,
            False,
            False,
            True,
            False,
            "sub",
        )

        assert "vertical-align:sub" in editor.toHtml().lower()
    finally:
        window._set_external_text_editor(None)
        window.close()


def test_external_rich_text_editor_preserves_selection_for_toolbar_action() -> None:
    window = ChemusonWindow()
    try:
        editor = QTextEdit()
        editor.setPlainText("O2")
        window._set_external_text_editor(editor)

        cursor = editor.textCursor()
        cursor.setPosition(1)
        cursor.movePosition(cursor.MoveOperation.End, cursor.MoveMode.KeepAnchor)
        editor.setTextCursor(cursor)
        window._external_text_cursor_state = (1, 2)
        window._external_text_selected_range = (1, 2)

        collapsed = editor.textCursor()
        collapsed.setPosition(2)
        editor.setTextCursor(collapsed)

        window.text_toolbar.action_sup.trigger()

        assert "vertical-align:super" in editor.toHtml().lower()
        assert window.text_toolbar.action_sup.shortcuts()
    finally:
        window._set_external_text_editor(None)
        window.close()


def test_ctrl_k_shortcut_triggers_full_clean_2d_action() -> None:
    window = ChemusonWindow()
    try:
        triggered: list[bool] = []
        window.action_clean_2d_full.triggered.connect(
            lambda _checked=False: triggered.append(True)
        )
        window.show()
        QApplication.processEvents()
        window.canvas.setFocus()

        QTest.keyClick(
            window.canvas.viewport(),
            Qt.Key.Key_K,
            Qt.KeyboardModifier.ControlModifier,
        )
        QApplication.processEvents()

        assert window.action_clean_2d_full.shortcut() == QKeySequence("Ctrl+K")
        assert triggered
    finally:
        window.close()


def test_ctrl_alt_k_shortcut_triggers_clean_2d_propose_action() -> None:
    window = ChemusonWindow()
    try:
        triggered: list[bool] = []
        window.action_clean_2d_propose.triggered.connect(
            lambda _checked=False: triggered.append(True)
        )
        window.show()
        QApplication.processEvents()
        window.canvas.setFocus()

        QTest.keyClick(
            window.canvas.viewport(),
            Qt.Key.Key_K,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.AltModifier,
        )
        QApplication.processEvents()

        assert window.action_clean_2d_propose.shortcut() == QKeySequence("Ctrl+Alt+K")
        assert triggered
    finally:
        window.close()
