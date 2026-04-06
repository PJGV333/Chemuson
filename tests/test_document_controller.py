"""Pruebas no-GUI para DocumentController."""

import os

from chemuson.gui.controllers.document_controller import (
    DocumentController,
    DocumentTabsContext,
    RecentFilesContext,
)


class _Settings:
    def __init__(self):
        self.store = {}

    def value(self, key, default=None):
        return self.store.get(key, default)

    def setValue(self, key, value):
        self.store[key] = value


class _UndoStack:
    def __init__(self, clean=True):
        self._clean = clean

    def isClean(self):
        return self._clean


class _Canvas:
    def __init__(self, clean=True):
        self.undo_stack = _UndoStack(clean=clean)


class _TabManager:
    def __init__(self, canvas):
        self.canvas = canvas
        self.file_paths = {canvas: None}
        self.tab_titles = {canvas: "Sin titulo"}

    def set_canvas_file_path(self, canvas, path):
        self.file_paths[canvas] = path
        if path:
            self.tab_titles[canvas] = os.path.basename(path)

    def update_tab_title(self, canvas):
        if canvas not in self.tab_titles:
            self.tab_titles[canvas] = "Sin titulo"


class _Window:
    def __init__(self):
        self._settings = _Settings()
        self._recent_files = []
        self.canvas = _Canvas(clean=True)
        self._tab_manager = _TabManager(self.canvas)
        self._current_file_path = None
        self.recent_menu_updated = False

    def _save_recent_files(self):
        DocumentController.save_recent_files(self.recent_context())

    def _update_recent_menu(self):
        self.recent_menu_updated = True

    def recent_context(self):
        return RecentFilesContext(
            settings=self._settings,
            recent_files=self._recent_files,
            persist_recent_files=self._save_recent_files,
            refresh_recent_menu=self._update_recent_menu,
        )

    def tabs_context(self):
        return DocumentTabsContext(
            tab_manager=self._tab_manager,
            current_canvas=self.canvas,
            current_file_path_setter=self._set_current_file_path,
        )

    def _set_current_file_path(self, path):
        self._current_file_path = path


def test_load_and_save_recent_files_roundtrip() -> None:
    window = _Window()
    window._settings.setValue("recent_files", ["a.cmsn", "b.cmsn"])
    loaded = DocumentController.load_recent_files(window)
    assert loaded == ["a.cmsn", "b.cmsn"]
    window._recent_files = loaded
    DocumentController.save_recent_files(window.recent_context())
    assert window._settings.value("recent_files") == ["a.cmsn", "b.cmsn"]


def test_set_canvas_file_path_updates_internal_state() -> None:
    window = _Window()
    filepath = os.path.abspath("demo.cmsn")
    DocumentController.set_canvas_file_path(window.tabs_context(), window.canvas, filepath)
    assert window._tab_manager.file_paths[window.canvas] == filepath
    assert window._tab_manager.tab_titles[window.canvas] == "demo.cmsn"
    assert window._current_file_path == filepath


def test_add_recent_file_moves_item_to_front_and_limits_size() -> None:
    window = _Window()
    window._recent_files = [os.path.abspath(f"{i}.cmsn") for i in range(12)]
    latest = os.path.abspath("5.cmsn")
    DocumentController.add_recent_file(window.recent_context(), latest)
    assert window._recent_files[0] == latest
    assert len(window._recent_files) == 10
    assert window.recent_menu_updated is True
