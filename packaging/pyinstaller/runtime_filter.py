"""Keep host-owned Linux desktop runtimes outside frozen PyInstaller archives."""

from __future__ import annotations

from pathlib import Path
from typing import MutableSequence

HOST_RUNTIME_LIBRARY_PREFIXES = (
    "libstdc++.so",
    "libfontconfig.so",
    "libfreetype.so",
)


def is_host_runtime_binary(entry: object) -> bool:
    """Return whether a PyInstaller binary entry shadows a host integration ABI."""
    if not isinstance(entry, (tuple, list)):
        return False
    return any(
        Path(str(part)).name.startswith(HOST_RUNTIME_LIBRARY_PREFIXES)
        for part in entry[:2]
    )


def exclude_host_runtime_binaries(entries: MutableSequence) -> tuple[str, ...]:
    """Remove host-owned C++/font runtimes from a mutable PyInstaller TOC."""
    removed: list[str] = []
    for entry in tuple(entries):
        if not is_host_runtime_binary(entry):
            continue
        entries.remove(entry)
        removed.append(Path(str(entry[0])).name)
    return tuple(removed)
