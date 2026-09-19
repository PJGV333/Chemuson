from __future__ import annotations

import json
import os

from chemuson.resilience.recovery import (
    archive_autosave,
    list_autosave_entries,
    read_autosave_metadata,
)


def _write_autosave(path, *, original_path="/tmp/example.cmsn", timestamp="2026-09-19T04:00:00") -> None:
    path.write_text(
        json.dumps(
            {
                "application": "Chemuson",
                "autosave_metadata": {
                    "original_path": original_path,
                    "timestamp": timestamp,
                },
            }
        ),
        encoding="utf-8",
    )


def test_recovery_policy_reads_lists_and_archives_autosaves(tmp_path) -> None:
    pending = tmp_path / "autosave"
    pending.mkdir()
    first = pending / "first.json"
    second = pending / "second.json"
    invalid = pending / "invalid.json"
    _write_autosave(first, original_path="/tmp/first.cmsn")
    _write_autosave(second, original_path="/tmp/second.cmsn")
    invalid.write_text(json.dumps({"application": "Other"}), encoding="utf-8")
    os.utime(first, (1, 1))
    os.utime(second, (2, 2))

    assert read_autosave_metadata(str(first)) == {
        "autosave_path": str(first),
        "original_path": "/tmp/first.cmsn",
        "timestamp": "2026-09-19T04:00:00",
    }
    assert read_autosave_metadata(str(invalid)) is None
    entries = list_autosave_entries(str(pending))
    assert [entry["filename"] for entry in entries] == ["second.json", "first.json"]

    archived = archive_autosave(str(second), str(pending))
    assert archived == str(pending / "old" / "second.json")
    assert not second.exists()
    assert (pending / "old" / "second.json").exists()


def test_recovery_policy_handles_missing_directory_and_missing_metadata(tmp_path) -> None:
    assert list_autosave_entries(str(tmp_path / "missing")) == []
    malformed = tmp_path / "malformed.json"
    malformed.write_text("not-json", encoding="utf-8")
    assert read_autosave_metadata(str(malformed)) is None
