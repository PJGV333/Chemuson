"""Behavioral and import-isolation tests for the autosave manager."""

from __future__ import annotations

import ast
import builtins
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Callable

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
    def __init__(
        self,
        interval_ms: int,
        single_shot: bool,
        callback: Callable[[], None],
    ) -> None:
        self.started = 0
        self.stopped = 0
        self.interval_ms = interval_ms
        self.single_shot = single_shot
        self.callback = callback

    def start(self) -> None:
        self.started += 1

    def stop(self) -> None:
        self.stopped += 1


def make_timer_factory(captured: list[FakeTimer]):
    def make_timer(
        interval_ms: int,
        single_shot: bool,
        callback: Callable[[], None],
    ) -> FakeTimer:
        timer = FakeTimer(interval_ms, single_shot, callback)
        captured.append(timer)
        return timer

    return make_timer


def make_manager(
    serializer: Callable[[FakeDocument], dict[str, object]],
    *,
    clean: bool = False,
    backup_limit: int = 5,
):
    timers: list[FakeTimer] = []
    document = FakeDocument()
    undo_stack = FakeUndoStack(clean)
    manager = AutosaveManager(
        document,
        undo_stack,
        serializer,
        make_timer_factory(timers),
        backup_limit=backup_limit,
    )
    return manager, document, undo_stack, timers


def test_write_autosave_creates_json_with_metadata(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    serialized_documents: list[FakeDocument] = []

    def serialize(document: FakeDocument) -> dict[str, object]:
        serialized_documents.append(document)
        return {"application": "Chemuson", "canvas": {"value": "preserved"}}

    manager, document, _, timers = make_manager(serialize)
    source_path = tmp_path / "sample.cmsn"
    manager.set_original_path(str(source_path))

    autosave_path = manager._write_autosave()

    assert autosave_path is not None
    assert os.path.isfile(autosave_path)
    assert serialized_documents == [document]
    assert [(timer.interval_ms, timer.single_shot) for timer in timers] == [
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
    manager, _, _, _ = make_manager(lambda _document: {"application": "Chemuson"}, backup_limit=2)

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
    manager, document, _, _ = make_manager(lambda value: calls.append(value) or {}, clean=True)

    assert manager._write_autosave() is None
    assert calls == []
    assert document is not None


def test_serializer_error_propagates_without_writing_success(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))

    def fail(_document: FakeDocument) -> dict[str, object]:
        raise RuntimeError("serializer failed")

    manager, _, _, _ = make_manager(fail)

    with pytest.raises(RuntimeError, match="serializer failed"):
        manager._write_autosave()
    assert manager._list_document_autosaves() == []


def test_start_is_idempotent_for_dirty_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {})

    manager.start()
    manager.start()

    assert timers[0].started == 1
    assert timers[1].started == 1


def test_start_only_starts_periodic_timer_for_clean_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {}, clean=True)

    manager.start()

    assert timers[0].started == 1
    assert timers[1].started == 0


def test_stop_and_restart_preserve_timer_lifecycle(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {})

    manager.start()
    manager.stop()
    manager.stop()
    manager.start()

    assert [timer.stopped for timer in timers] == [2, 2]
    assert [timer.started for timer in timers] == [2, 2]


def test_restart_debounce_is_inactive_before_start(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {})

    manager.restart_debounce()

    assert timers[1].started == 0


def test_restart_debounce_starts_for_running_dirty_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {})
    manager.start()

    manager.restart_debounce()

    assert timers[1].started == 2


def test_restart_debounce_stops_for_running_clean_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, undo_stack, timers = make_manager(lambda _document: {})
    manager.start()
    undo_stack.clean = True

    manager.restart_debounce()

    assert timers[1].started == 1
    assert timers[1].stopped == 1


def test_cancel_debounce_stops_timer(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, timers = make_manager(lambda _document: {})

    manager.cancel_debounce()

    assert timers[1].stopped == 1


def test_periodic_timer_callback_writes_dirty_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    serialized: list[FakeDocument] = []
    manager, document, _, timers = make_manager(lambda value: serialized.append(value) or {"data": "kept"})

    timers[0].callback()

    assert serialized == [document]
    assert len(manager._list_document_autosaves()) == 1


def test_periodic_timer_callback_skips_clean_document(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    serialized: list[FakeDocument] = []
    manager, _, _, timers = make_manager(lambda value: serialized.append(value) or {}, clean=True)

    timers[0].callback()

    assert serialized == []
    assert manager._list_document_autosaves() == []


def test_cleanup_after_save_replaces_snapshots_with_forced_metadata_snapshot(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, _ = make_manager(lambda _document: {"application": "Chemuson"}, clean=True)
    manager.set_original_path(str(tmp_path / "saved.cmsn"))
    for _ in range(2):
        assert manager._write_autosave(force=True) is not None
        time.sleep(0.01)

    manager.cleanup_after_save()

    autosaves = manager._list_document_autosaves()
    assert len(autosaves) == 1
    with open(autosaves[0], encoding="utf-8") as snapshot:
        metadata = json.load(snapshot)["autosave_metadata"]
    assert metadata["original_path"] == os.path.abspath(str(tmp_path / "saved.cmsn"))
    assert metadata["doc_id"] == os.path.abspath(str(tmp_path / "saved.cmsn"))
    assert metadata["doc_hash"] == manager._doc_hash()
    assert metadata["timestamp"]


def test_document_hash_is_stable_and_hexadecimal(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    path = str(tmp_path / "stable.cmsn")
    first, _, _, _ = make_manager(lambda _document: {})
    second, _, _, _ = make_manager(lambda _document: {})
    first.set_original_path(path)
    second.set_original_path(path)

    assert first._doc_hash() == second._doc_hash()
    assert re.fullmatch(r"[0-9a-f]{16}", first._doc_hash())


def test_autosave_filename_uses_stable_recovery_format(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, _ = make_manager(lambda _document: {})

    assert re.fullmatch(r"[0-9a-f]{16}_[0-9]{8}_[0-9]{6}_[0-9]{6}\.json", manager._autosave_name())


def test_metadata_preserves_recovery_fields_and_serializer_payload(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, _ = make_manager(lambda _document: {"application": "Chemuson", "canvas": {"value": 1}})
    manager.set_original_path(str(tmp_path / "metadata.cmsn"))

    autosave_path = manager._write_autosave()

    assert autosave_path is not None
    with open(autosave_path, encoding="utf-8") as snapshot:
        payload = json.load(snapshot)
    assert payload["application"] == "Chemuson"
    assert payload["canvas"] == {"value": 1}
    assert set(payload["autosave_metadata"]) == {"original_path", "timestamp", "doc_id", "doc_hash"}
    assert payload["autosave_metadata"]["original_path"] == os.path.abspath(str(tmp_path / "metadata.cmsn"))
    assert payload["autosave_metadata"]["doc_id"] == os.path.abspath(str(tmp_path / "metadata.cmsn"))
    assert payload["autosave_metadata"]["doc_hash"] == manager._doc_hash()
    assert payload["autosave_metadata"]["timestamp"]


def test_write_error_propagates_without_creating_snapshot(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    manager, _, _, _ = make_manager(lambda _document: {"application": "Chemuson"})

    def fail_open(*_args, **_kwargs):
        raise OSError("write failed")

    monkeypatch.setattr(builtins, "open", fail_open)

    with pytest.raises(OSError, match="write failed"):
        manager._write_autosave()
    assert manager._list_document_autosaves() == []


def is_blocked_module_name(name: str, blocked: tuple[str, ...]) -> bool:
    return any(name == prefix or name.startswith(prefix + ".") for prefix in blocked)


def test_autosave_module_ast_has_no_gui_chemio_qt_or_rdkit_imports() -> None:
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
        if is_blocked_module_name(name, ("PyQt6", "chemuson.gui", "chemuson.chemio", "rdkit"))
    ]


def test_blocked_module_matching_catches_rdkit_submodules() -> None:
    blocked = ("PyQt6", "chemuson.gui", "chemuson.chemio", "rdkit")

    assert is_blocked_module_name("rdkit.Chem", blocked)
    assert not is_blocked_module_name("rdkittools", blocked)


def test_autosave_import_isolated_in_subprocess() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import chemuson.utils.autosave; "
            "blocked = ('PyQt6', 'chemuson.gui', 'chemuson.chemio', 'rdkit'); "
            "unexpected = [name for name in sys.modules if any(name == prefix or name.startswith(prefix + '.') for prefix in blocked)]; "
            "print('unexpected modules:', unexpected); raise SystemExit(bool(unexpected))",
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
    assert result.returncode == 0, result.stdout + result.stderr
