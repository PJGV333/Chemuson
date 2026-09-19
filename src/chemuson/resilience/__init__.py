"""Resilience services for crash logging and document recovery."""

from .autosave import AutosaveManager
from .crash_reporter import install, write_crash_log

__all__ = ["AutosaveManager", "install", "write_crash_log"]
