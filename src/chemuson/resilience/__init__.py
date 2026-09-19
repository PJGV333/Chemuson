"""Resilience services for crash logging and document recovery."""

from .autosave import AutosaveManager
from .crash_reporter import install, write_crash_log
from .recovery import archive_autosave, list_autosave_entries, read_autosave_metadata

__all__ = [
    "AutosaveManager",
    "archive_autosave",
    "install",
    "list_autosave_entries",
    "read_autosave_metadata",
    "write_crash_log",
]
