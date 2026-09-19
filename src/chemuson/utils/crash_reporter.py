"""Compatibility shim for the historical crash reporter path."""

from chemuson.resilience.crash_reporter import install, write_crash_log

__all__ = ["install", "write_crash_log"]
