"""Pruebas no-GUI para DocumentController."""

import os

from chemuson.gui.controllers.document_controller import DocumentController


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


class _AutosaveManager:
    def __init__(self):
        self.original_path = None

    def set_original_path(self, path):
        self.original_path = path


class _Tabs:
    def __init__(self, canvas):
        self._canvas = canvas
        self._text = "Documento"
        self._tooltip = ""

    def indexOf(self, canvas):
        return 0 if canvas is self._canvas else -1

    def setTabText(self, idx, text):
        self._text = text

    def setTabToolTip(self, idx, text):
        self._tooltip = text


class _Window:
    def __init__(self):
        self._settings = _Settings()
        self._recent_files = []
        self.canvas = _Canvas(clean=True)
        self._canvas_file_paths = {self.canvas: None}
        self._canvas_tab_titles = {self.canvas: "Sin título"}
        self._canvas_autosave_managers = {self.canvas: _AutosaveManager()}
        self._current_file_path = None
        self.tabs = _Tabs(self.canvas)
        self.recent_menu_updated = False

    def _tabs_alive(self):
        return True

    def _save_recent_files(self):
        DocumentController.save_recent_files(self)

    def _update_recent_menu(self):
        self.recent_menu_updated = True


def test_load_and_save_recent_files_roundtrip() -> None:
    window = _Window()
    window._settings.setValue("recent_files", ["a.cmsn", "b.cmsn"])
    loaded = DocumentController.load_recent_files(window)
    assert loaded == ["a.cmsn", "b.cmsn"]
    window._recent_files = loaded
    DocumentController.save_recent_files(window)
    assert window._settings.value("recent_files") == ["a.cmsn", "b.cmsn"]


def test_set_canvas_file_path_updates_internal_state() -> None:
    window = _Window()
    filepath = os.path.abspath("demo.cmsn")
    DocumentController.set_canvas_file_path(window, window.canvas, filepath)
    assert window._canvas_file_paths[window.canvas] == filepath
    assert window._canvas_tab_titles[window.canvas] == "demo.cmsn"
    assert window._current_file_path == filepath
    assert window._canvas_autosave_managers[window.canvas].original_path == filepath


def test_add_recent_file_moves_item_to_front_and_limits_size() -> None:
    window = _Window()
    window._recent_files = [os.path.abspath(f"{i}.cmsn") for i in range(12)]
    latest = os.path.abspath("5.cmsn")
    DocumentController.add_recent_file(window, latest)
    assert window._recent_files[0] == latest
    assert len(window._recent_files) == 10
    assert window.recent_menu_updated is True
