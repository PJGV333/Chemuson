"""Pruebas para el script de sincronización de versión de release."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load_module():
    script_path = (
        Path(__file__).resolve().parent.parent
        / "packaging"
        / "release"
        / "set_version.py"
    )
    spec = importlib.util.spec_from_file_location("set_version", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_write_version_updates_version_file(tmp_path) -> None:
    module = _load_module()
    version_file = tmp_path / "_version.py"
    version_file.write_text('__version__ = "0.2.1"\n', encoding="utf-8")

    normalized = module.write_version(version_file, "0.2.2-beta.1")

    assert normalized == "0.2.2-beta.1"
    assert version_file.read_text(encoding="utf-8") == '__version__ = "0.2.2-beta.1"\n'


def test_validate_version_rejects_invalid_semver() -> None:
    module = _load_module()

    with pytest.raises(ValueError):
        module.validate_version("hotfix-now")
