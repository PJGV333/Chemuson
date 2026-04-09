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


def test_update_flatpak_metainfo_syncs_version_license_and_release_date(tmp_path) -> None:
    module = _load_module()
    metainfo_file = tmp_path / "io.github.PJGV333.Chemuson.metainfo.xml"
    metainfo_file.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<component type="desktop-application">
  <id>io.github.PJGV333.Chemuson</id>
  <project_license>MIT</project_license>
  <releases>
    <release version="0.2.1" date="2026-02-26"/>
  </releases>
</component>
""",
        encoding="utf-8",
    )

    normalized, normalized_date = module.update_flatpak_metainfo(
        metainfo_file=metainfo_file,
        version="0.2.3-beta.2",
        release_date="2026-04-09",
    )

    text = metainfo_file.read_text(encoding="utf-8")
    assert normalized == "0.2.3-beta.2"
    assert normalized_date == "2026-04-09"
    assert "<project_license>AGPL-3.0-only</project_license>" in text
    assert '<release version="0.2.3-beta.2" date="2026-04-09" />' in text
    assert '<release version="0.2.1" date="2026-02-26" />' in text
    assert text.index('version="0.2.3-beta.2"') < text.index('version="0.2.1"')


def test_validate_version_rejects_invalid_semver() -> None:
    module = _load_module()

    with pytest.raises(ValueError):
        module.validate_version("hotfix-now")


def test_validate_release_date_rejects_invalid_iso_date() -> None:
    module = _load_module()

    with pytest.raises(ValueError):
        module.validate_release_date("09-04-2026")
