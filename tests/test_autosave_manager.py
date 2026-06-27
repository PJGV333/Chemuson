"""Pruebas unitarias para el gestor de autosaves."""

import json
import os
import time

import pytest
from PyQt6.QtGui import QUndoCommand
from PyQt6.QtWidgets import QApplication


from chemuson.gui.canvas import ChemusonCanvas
from chemuson.utils.autosave import AutosaveManager


class _NoOpCommand(QUndoCommand):
    """Comando mínimo para ensuciar el undo stack en pruebas."""

    def redo(self) -> None:
        return

    def undo(self) -> None:
        return


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    return QApplication.instance() or QApplication([])


def _make_dirty_canvas() -> ChemusonCanvas:
    canvas = ChemusonCanvas()
    canvas.undo_stack.push(_NoOpCommand("dirty"))
    return canvas


def test_write_autosave_creates_json_with_metadata(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    canvas = _make_dirty_canvas()
    manager = AutosaveManager(canvas, canvas, canvas.undo_stack)
    source_path = tmp_path / "sample.cmsn"
    manager.set_original_path(str(source_path))

    autosave_path = manager._write_autosave()

    assert autosave_path is not None
    assert os.path.isfile(autosave_path)
    with open(autosave_path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    metadata = payload.get("autosave_metadata", {})
    assert metadata.get("original_path") == os.path.abspath(str(source_path))
    assert metadata.get("timestamp")
    assert payload.get("application") == "Chemuson"


def test_rotation_removes_oldest_autosaves(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    canvas = _make_dirty_canvas()
    manager = AutosaveManager(canvas, canvas, canvas.undo_stack, backup_limit=2)

    created_paths = []
    for _ in range(4):
        autosave_path = manager._write_autosave()
        assert autosave_path is not None
        created_paths.append(autosave_path)
        time.sleep(0.01)

    remaining = manager._list_document_autosaves()
    assert len(remaining) == 2
    assert os.path.exists(created_paths[-1])
    assert os.path.exists(created_paths[-2])
    assert not os.path.exists(created_paths[0])
    assert not os.path.exists(created_paths[1])
