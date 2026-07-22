"""Behavioral and import-isolation tests for the autosave manager."""

from __future__ import annotations

import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from chemuson.utils.autosave import AutosaveManager


class FakeDocument:
    pass


class FakeUndoStack:
    def __init__(self, clean: bool = False) -> None:
        self.clean = clean

    def isClean(self) -> bool:
        return self.clean


class FakeTimer:
    def __init__(self) -> None:
        self.started = 0
        self.stopped = 0

    def start(self) -> None:
        self.started += 1

    def stop(self) -> None:
        self.stopped += 1


def make_timer_factory(captured: list[tuple[int, bool, object]]) -> FakeTimer:
    def make_timer(interval_ms: int, single_shot: bool, callback: object) -> FakeTimer:
        captured.append((interval_ms, single_shot, callback))
        return FakeTimer()

    return make_timer


def make_manager(serializer, *, clean: bool = False, backup_limit: int = 5):
    timer_specs: list[tuple[int, bool, object]] = []
    document = FakeDocument()
    manager = AutosaveManager(
        document,
        FakeUndoStack(clean),
        serializer,
        make_timer_factory(timer_specs),
        backup_limit=backup_limit,
    )
    return manager, document, timer_specs


def test_write_autosave_creates_json_with_metadata(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    serialized_documents: list[FakeDocument] = []

    def serialize(document: FakeDocument) -> dict[str, object]:
        serialized_documents.append(document)
        return {"application": "Chemuson", "canvas": {"value": "preserved"}}

    manager, document, timer_specs = make_manager(serialize)
    source_path = tmp_path / "sample.cmsn"
    manager.set_original_path(str(source_path))

    autosave_path = manager._write_autosave()

    assert autosave_path is not None
    assert os.path.isfile(autosave_path)
    assert serialized_documents == [document]
    assert [(interval, single_shot) for interval, single_shot, _ in timer_specs] == [
        (AutosaveManager.AUTOSAVE_INTERVAL_MS, False),
        (AutosaveManager.DEBOUNCE_INTERVAL_MS, True),
    ]
    with open(autosave_path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    metadata = payload["autosave_metadata"]
    assert metadata["original_path"] == os.path.abspath(str(source_path))
    assert metadata["timestamp"]
    assert payload["application"] == "Chemuson"
    assert payload["canvas"] == {"value": "preserved"}


def test_rotation_removes_oldest_autosaves(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _ = make_manager(lambda _document: {"application": "Chemuson"}, backup_limit=2)

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


def test_clean_document_does_not_invoke_serializer(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    calls: list[FakeDocument] = []
    manager, document, _ = make_manager(lambda value: calls.append(value) or {}, clean=True)

    assert manager._write_autosave() is None
    assert calls == []
    assert document is not None


def test_serializer_error_propagates_without_writing_success(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))

    def fail(_document: FakeDocument) -> dict[str, object]:
        raise RuntimeError("serializer failed")

    manager, _, _ = make_manager(fail)

    with pytest.raises(RuntimeError, match="serializer failed"):
        manager._write_autosave()
    assert manager._list_document_autosaves() == []


def test_autosave_module_ast_has_no_gui_chemio_or_qt_imports() -> None:
    path = Path(__file__).resolve().parents[1] / "src/chemuson/utils/autosave.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imported_modules = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.append(node.module or "")
    assert not [
        name
        for name in imported_modules
        if name == "PyQt6"
        or name.startswith("PyQt6.")
        or name == "chemuson.gui"
        or name.startswith("chemuson.gui.")
        or name == "chemuson.chemio"
        or name.startswith("chemuson.chemio.")
    ]


def test_autosave_import_isolated_in_subprocess() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import chemuson.utils.autosave; "
            "blocked = ('PyQt6', 'chemuson.gui', 'chemuson.chemio'); "
            "raise SystemExit(bool([name for name in sys.modules if any(name == prefix or name.startswith(prefix + '.') for prefix in blocked)]))",
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
        env={
            **os.environ,
            "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src"),
        },
    )
    assert result.returncode == 0, result.stderr
