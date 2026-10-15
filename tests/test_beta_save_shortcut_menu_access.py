"""Interacción de guardado, cierre y modos de menú bloqueantes de la beta."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from PyQt6.QtCore import Qt, QSettings
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox

from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.main_window import ChemusonWindow
from chemuson.platform.settings import load_ui_preferences


@pytest.fixture
def window(tmp_path, monkeypatch):
    settings_path = tmp_path / "Chemuson.ini"

    def settings_factory():
        return QSettings(str(settings_path), QSettings.Format.IniFormat)

    monkeypatch.setattr("chemuson.platform.settings.application_settings", settings_factory)
    monkeypatch.setattr("chemuson.gui.shell.assembly.application_settings", settings_factory)
    settings = settings_factory()
    settings.clear()
    settings.setValue("ui/onboarding/completed", True)
    settings.sync()

    instance = ChemusonWindow()
    instance.show()
    QApplication.processEvents()
    QApplication.setActiveWindow(instance)
    instance.canvas.setFocus()
    QApplication.processEvents()
    yield instance

    for index in range(instance.tabs.count()):
        canvas = instance._canvas_from_tab_index(index)
        if canvas is not None:
            canvas.undo_stack.setClean()
    instance.close()
    QApplication.processEvents()


def _save_path(window: ChemusonWindow, canvas, path: Path) -> None:
    PersistenceManager.save_to_file(str(path), canvas)
    window._tab_manager.set_canvas_file_path(canvas, str(path))


def _add_atom_and_mark_dirty(canvas, element: str) -> None:
    canvas.model.add_atom(element, 12.0, 18.0)
    canvas._rebuild_items_from_model()
    canvas.undo_stack.resetClean()
    assert not canvas.undo_stack.isClean()


def _saved_elements(path: Path) -> list[str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return [atom["element"] for atom in payload["model"]["atoms"]]


def _focus(window: ChemusonWindow, canvas) -> None:
    window.tabs.setCurrentWidget(canvas)
    QApplication.setActiveWindow(window)
    canvas.setFocus()
    QApplication.processEvents()


def _close_button(window: ChemusonWindow, canvas):
    index = window.tabs.indexOf(canvas)
    assert index >= 0
    side = window.app_bar.tab_bar._side_button(index)
    assert side is not None
    assert side.close_button.isVisible()
    return side.close_button


@pytest.mark.parametrize("classic_menu", [False, True], ids=["compact", "classic"])
def test_ctrl_s_saves_only_active_existing_document_through_original_action(
    window, tmp_path, monkeypatch, classic_menu
) -> None:
    if classic_menu:
        window.action_show_menu_bar.trigger()
    assert window.menuBar().isVisible() is classic_menu
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *_args: pytest.fail("existing CMSN path prompted for a route")),
    )
    active = window.canvas
    active_path = tmp_path / "active.cmsn"
    _save_path(window, active, active_path)
    _add_atom_and_mark_dirty(active, "C")

    window._on_file_new()
    inactive = window.canvas
    inactive_path = tmp_path / "inactive.cmsn"
    _save_path(window, inactive, inactive_path)
    _add_atom_and_mark_dirty(inactive, "N")
    _focus(window, active)

    file_menu = next(
        action.menu()
        for action in window.menuBar().actions()
        if action.menu() is not None and action.menu().title() == "Archivo"
    )
    assert window.action_save in file_menu.actions()
    assert window.action_save in window.actions()
    assert [shortcut.toString() for shortcut in window.action_save.shortcuts()] == [
        "Ctrl+S",
        "Ctrl+G",
    ]
    triggered: list[bool] = []
    window.action_save.triggered.connect(lambda _checked=False: triggered.append(True))
    controller_calls: list[object] = []
    save_file = window._file_controller.save_file
    monkeypatch.setattr(
        window._file_controller,
        "save_file",
        lambda context: (controller_calls.append(context.canvas), save_file(context)),
    )

    QTest.keyClick(active, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()

    assert triggered == [True]
    assert controller_calls == [active]
    assert _saved_elements(active_path) == ["C"]
    assert _saved_elements(inactive_path) == []
    assert active.undo_stack.isClean()
    assert not inactive.undo_stack.isClean()
    assert window.canvas is active


def test_ctrl_g_remains_an_alias_for_the_same_save_action(window, tmp_path) -> None:
    canvas = window.canvas
    path = tmp_path / "ctrl-g.cmsn"
    _save_path(window, canvas, path)
    _add_atom_and_mark_dirty(canvas, "O")
    triggered: list[bool] = []
    window.action_save.triggered.connect(lambda _checked=False: triggered.append(True))

    QTest.keyClick(canvas, Qt.Key.Key_G, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()

    assert triggered == [True]
    assert _saved_elements(path) == ["O"]
    assert canvas.undo_stack.isClean()


def test_menu_save_uses_same_action_and_writes_cmsn(window, tmp_path) -> None:
    canvas = window.canvas
    path = tmp_path / "menu.cmsn"
    _save_path(window, canvas, path)
    _add_atom_and_mark_dirty(canvas, "F")

    window.action_show_menu_bar.trigger()
    QApplication.processEvents()
    assert window.menuBar().isVisible()
    file_action = next(
        action
        for action in window.menuBar().actions()
        if action.menu() is not None and action.menu().title() == "Archivo"
    )
    file_menu = file_action.menu()
    QTest.mouseClick(
        window.menuBar(),
        Qt.MouseButton.LeftButton,
        pos=window.menuBar().actionGeometry(file_action).center(),
    )
    QApplication.processEvents()
    assert file_menu.isVisible()
    QTest.mouseClick(
        file_menu,
        Qt.MouseButton.LeftButton,
        pos=file_menu.actionGeometry(window.action_save).center(),
    )
    QApplication.processEvents()

    assert _saved_elements(path) == ["F"]
    assert canvas.undo_stack.isClean()


def test_ctrl_s_prompts_for_a_new_native_path_and_saves_content(
    window, tmp_path, monkeypatch
) -> None:
    canvas = window.canvas
    _add_atom_and_mark_dirty(canvas, "Cl")
    destination = tmp_path / "new-document.cmsn"
    calls: list[tuple] = []
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *args: (calls.append(args) or (str(destination), "Documento ChemUSON (*.cmsn)"))),
    )

    QTest.keyClick(canvas, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()

    assert len(calls) == 1
    assert destination.is_file()
    assert _saved_elements(destination) == ["Cl"]
    assert window._tab_manager.file_path_for(canvas) == str(destination.resolve())
    assert canvas.undo_stack.isClean()


def test_cancelled_new_document_save_keeps_dirty_content_and_path(
    window, monkeypatch
) -> None:
    canvas = window.canvas
    _add_atom_and_mark_dirty(canvas, "Br")
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *_args: ("", "")),
    )

    QTest.keyClick(canvas, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()

    assert len(canvas.model.atoms) == 1
    assert not canvas.undo_stack.isClean()
    assert window._tab_manager.file_path_for(canvas) is None


def test_failed_existing_file_save_keeps_original_file_and_dirty_state(
    window, tmp_path, monkeypatch
) -> None:
    canvas = window.canvas
    path = tmp_path / "protected.cmsn"
    _save_path(window, canvas, path)
    original_bytes = path.read_bytes()
    _add_atom_and_mark_dirty(canvas, "S")

    def fail_save(*_args, **_kwargs):
        raise OSError("simulated write failure")

    monkeypatch.setattr(PersistenceManager, "save_to_file", fail_save)
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *_args: None))
    QTest.keyClick(canvas, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()

    assert path.read_bytes() == original_bytes
    assert len(canvas.model.atoms) == 1
    assert window._tab_manager.file_path_for(canvas) == str(path.resolve())
    assert not canvas.undo_stack.isClean()


def _dirty_nonactive_tab(window, tmp_path):
    target = window.canvas
    path = tmp_path / "nonactive.cmsn"
    _save_path(window, target, path)
    _add_atom_and_mark_dirty(target, "P")
    window._on_file_new()
    remaining = window.canvas
    return target, remaining, path


def test_close_button_save_writes_nonactive_target_before_closing(
    window, tmp_path, monkeypatch
) -> None:
    target, remaining, path = _dirty_nonactive_tab(window, tmp_path)
    replies: list[str] = []
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(
            lambda *_args: (replies.append("asked") or QMessageBox.StandardButton.Save)
        ),
    )

    QTest.mouseClick(_close_button(window, target), Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert replies == ["asked"]
    assert _saved_elements(path) == ["P"]
    assert window.tabs.indexOf(target) == -1
    assert window.tabs.count() == 1
    assert window.canvas is remaining
    assert remaining.undo_stack.isClean()


@pytest.mark.parametrize(
    ("reply", "closed"),
    [
        pytest.param(QMessageBox.StandardButton.Discard, True, id="discard"),
        pytest.param(QMessageBox.StandardButton.Cancel, False, id="cancel"),
    ],
)
def test_close_button_discard_or_cancel_is_explicit_and_only_affects_target(
    window, tmp_path, monkeypatch, reply, closed
) -> None:
    target, remaining, path = _dirty_nonactive_tab(window, tmp_path)
    questions: list[str] = []
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *_args: (questions.append("asked") or reply)),
    )

    QTest.mouseClick(_close_button(window, target), Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert questions == ["asked"]
    if closed:
        assert window.tabs.indexOf(target) == -1
        assert window.tabs.count() == 1
    else:
        assert window.tabs.indexOf(target) >= 0
        assert window.tabs.count() == 2
        assert len(target.model.atoms) == 1
        assert not target.undo_stack.isClean()
    assert window.canvas is remaining
    assert _saved_elements(path) == []


def test_clean_tab_closes_without_discard_prompt(window, monkeypatch) -> None:
    clean = window.canvas
    window._on_file_new()
    active = window.canvas
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *_args: pytest.fail("clean tab prompted for discard")),
    )

    QTest.mouseClick(_close_button(window, clean), Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert window.tabs.count() == 1
    assert window.tabs.indexOf(clean) == -1
    assert window.canvas is active


def test_close_save_cancelled_new_route_keeps_dirty_tab_open(
    window, tmp_path, monkeypatch
) -> None:
    target = window.canvas
    _add_atom_and_mark_dirty(target, "P")
    window._on_file_new()
    remaining = window.canvas
    replies: list[str] = []
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(
            lambda *_args: (replies.append("asked") or QMessageBox.StandardButton.Save)
        ),
    )
    dialog_calls: list[bool] = []
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *_args: (dialog_calls.append(True) or ("", ""))),
    )

    QTest.mouseClick(_close_button(window, target), Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert replies == ["asked"]
    assert dialog_calls == [True]
    assert window.tabs.indexOf(target) >= 0
    assert window.tabs.count() == 2
    assert window.canvas is remaining
    assert len(target.model.atoms) == 1
    assert not target.undo_stack.isClean()


def test_close_save_write_failure_keeps_dirty_tab_open(
    window, tmp_path, monkeypatch
) -> None:
    target, remaining, path = _dirty_nonactive_tab(window, tmp_path)
    original_bytes = path.read_bytes()
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *_args: QMessageBox.StandardButton.Save),
    )
    monkeypatch.setattr(
        PersistenceManager,
        "save_to_file",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *_args: None))

    QTest.mouseClick(_close_button(window, target), Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert path.read_bytes() == original_bytes
    assert window.tabs.indexOf(target) >= 0
    assert window.tabs.count() == 2
    assert window.canvas is remaining
    assert len(target.model.atoms) == 1
    assert not target.undo_stack.isClean()


def test_menu_mode_toggle_persists_and_preserves_layout_across_themes(
    window, tmp_path
) -> None:
    assert not window.action_show_menu_bar.isChecked()
    assert not window.menuBar().isVisible()
    assert window.action_show_menu_bar.text() == "Mostrar barra de menús"
    assert window.action_show_menu_bar in window.view_menu.actions()
    aggregated_view = next(
        action.menu()
        for action in window._menu_aggregator.actions()
        if action.menu() is not None and action.menu().title() == "Ver"
    )
    assert aggregated_view is window.view_menu

    window.action_show_menu_bar.trigger()
    QApplication.processEvents()
    assert window.action_show_menu_bar.isChecked()
    assert window.menuBar().isVisible()
    prefs = load_ui_preferences(window._settings)
    assert prefs.show_menu_bar is True

    expected_titles = ("Archivo", "Editar", "Ver", "Estructura", "Reacción", "Ayuda")
    original_menus = [action.menu() for action in window.menuBar().actions()]
    original_save = window.action_save
    original_canvas = window.tabs
    original_palette = window.command_palette
    for theme in ("light", "dark"):
        window.current_theme = theme
        window._apply_theme()
        for width, height in ((1440, 900), (980, 600)):
            window.resize(width, height)
            QApplication.processEvents()
            assert window.menuBar().isVisible()
            assert [menu.title() for menu in original_menus] == list(expected_titles)
            assert [action.menu() for action in window.menuBar().actions()] == original_menus
            assert window.action_save is original_save
            assert window.tabs is original_canvas
            assert window.command_palette is original_palette

    window.action_show_menu_bar.trigger()
    QApplication.processEvents()
    assert not window.menuBar().isVisible()
    assert window.app_bar.menu_button.isVisible()
    assert load_ui_preferences(window._settings).show_menu_bar is False

    # La preferencia guardada se aplica en una nueva ventana.
    window.action_show_menu_bar.trigger()
    QApplication.processEvents()
    second = ChemusonWindow()
    try:
        second.show()
        QApplication.processEvents()
        assert second.action_show_menu_bar.isChecked()
        assert second.menuBar().isVisible()
    finally:
        for index in range(second.tabs.count()):
            canvas = second._canvas_from_tab_index(index)
            if canvas is not None:
                canvas.undo_stack.setClean()
        second.close()
        QApplication.processEvents()


def test_save_then_theme_preference_save_does_not_reset_menu_mode(window) -> None:
    window.action_show_menu_bar.trigger()
    window._apply_preferences({"theme": "dark"})
    assert window.menuBar().isVisible()
    assert load_ui_preferences(window._settings).show_menu_bar is True
