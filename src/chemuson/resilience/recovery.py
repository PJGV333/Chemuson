"""GUI-free filesystem policy for autosave recovery entries."""

from datetime import datetime
import json
import os
from typing import Optional


def read_autosave_metadata(filepath: str) -> Optional[dict]:
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            payload = json.load(f)
    except Exception:
        return None
    if payload.get("application") != "Chemuson":
        return None
    metadata = payload.get("autosave_metadata")
    if not isinstance(metadata, dict):
        metadata = {}
    raw_path = metadata.get("original_path")
    raw_timestamp = metadata.get("timestamp")
    return {
        "autosave_path": filepath,
        "original_path": str(raw_path) if raw_path else None,
        "timestamp": str(raw_timestamp) if raw_timestamp else "Desconocida",
    }


def list_autosave_entries(directory: str) -> list[dict]:
    if not os.path.isdir(directory):
        return []
    entries: list[dict] = []
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".json"):
            continue
        filepath = os.path.join(directory, name)
        if not os.path.isfile(filepath):
            continue
        metadata = read_autosave_metadata(filepath)
        if metadata is None:
            continue
        metadata["filename"] = name
        entries.append(metadata)
    entries.sort(key=lambda entry: os.path.getmtime(entry["autosave_path"]), reverse=True)
    return entries


def archive_autosave(path: str, autosave_dir: str) -> str:
    old_dir = os.path.join(autosave_dir, "old")
    os.makedirs(old_dir, exist_ok=True)
    basename = os.path.basename(path)
    target = os.path.join(old_dir, basename)
    if os.path.exists(target):
        root, ext = os.path.splitext(basename)
        target = os.path.join(old_dir, f"{root}_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}{ext}")
    os.replace(path, target)
    return target


__all__ = [
    "archive_autosave",
    "list_autosave_entries",
    "read_autosave_metadata",
]
