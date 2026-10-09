"""Regression tests for visible ChemUSON branding and stable technical IDs."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _text(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def test_about_description_and_main_window_title_use_approved_copy() -> None:
    about_source = _text("src/chemuson/gui/main_window.py")
    shell_source = _text("src/chemuson/gui/shell/assembly.py")

    assert "ChemUSON es un editor molecular libre y de código abierto para " in about_source
    assert "crear, editar, visualizar y analizar estructuras químicas, diagramas " in about_source
    assert "y anotaciones científicas." in about_source
    assert "inspirado en ChemDoodle" not in about_source
    assert (
        'self.setWindowTitle(f"ChemUSON {self._app_version} — Editor Molecular Libre")'
        in shell_source
    )


def test_linux_launchers_and_appstream_present_official_name() -> None:
    flatpak_desktop = _text("packaging/flatpak/io.github.PJGV333.Chemuson.desktop")
    appimage_desktop = _text(
        "packaging/linux/appimage/io.github.PJGV333.Chemuson.desktop.in"
    )
    metainfo_path = ROOT / "packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml"
    metainfo = ET.parse(metainfo_path).getroot()

    assert "Name=ChemUSON\n" in flatpak_desktop
    assert "Comment=ChemUSON —" in flatpak_desktop
    assert "Name=ChemUSON\n" in appimage_desktop
    assert "Comment=ChemUSON —" in appimage_desktop
    assert metainfo.findtext("name") == "ChemUSON"
    assert metainfo.findtext("id") == "io.github.PJGV333.Chemuson"


def test_windows_installer_display_name_changes_without_changing_upgrade_identity() -> None:
    installer = _text("packaging/windows/Chemuson.iss")

    assert '#define MyAppName "ChemUSON"' in installer
    assert '#define MyAppPublisher "ChemUSON"' in installer
    assert "DefaultGroupName=ChemUSON" in installer
    assert 'AppId={{E2D4477C-AE35-4C8E-9F7E-8C2E4DBE69A7}' in installer
    assert '#define MyAppExeName "Chemuson.exe"' in installer
    assert r"DefaultDirName={autopf}\Chemuson" in installer
    assert "OutputBaseFilename=Chemuson-v{#MyAppVersion}-windows-x86_64-setup" in installer


def test_package_import_settings_update_and_project_ids_remain_compatible() -> None:
    version = _text("src/chemuson/_version.py")
    pyproject = _text("pyproject.toml")
    settings = _text("src/chemuson/platform/settings.py")
    persistence = _text("src/chemuson/chemio/persistence.py")
    update_controller = _text("src/chemuson/gui/controllers/update_controller.py")
    bootstrap = _text("src/chemuson/app/bootstrap.py")
    flatpak_manifest = _text("packaging/flatpak/io.github.PJGV333.Chemuson.yml")
    preview_workflow = _text(".github/workflows/build-preview.yml")
    appimage_validator = _text("packaging/release/validate_appimage.py")
    pages_index = _text("packaging/release/generate_flatpak_pages_index.py")

    assert '__version__ = "0.3.0-beta.1"' in version
    assert re.search(r"(?m)^name\s*=\s*['\"]chemuson['\"]", pyproject)
    assert 'QSettings("Chemuson", "Chemuson")' in settings
    assert 'app.setApplicationName("Chemuson")' in bootstrap
    assert 'app.setApplicationDisplayName("ChemUSON")' in bootstrap
    assert '"application": "Chemuson"' in persistence
    assert 'app-id: io.github.PJGV333.Chemuson' in flatpak_manifest
    assert 'GitHubReleasesProvider("PJGV333", "Chemuson"' in update_controller
    assert "chemuson-preview-windows-portable" in preview_workflow
    assert "chemuson-preview-linux-flatpak" in preview_workflow
    assert 'entry.get("Name") != "ChemUSON"' in appimage_validator
    assert "<title>ChemUSON Flatpak</title>" in pages_index
    assert 'default="Chemuson"' in pages_index


def test_flatpak_remote_display_title_is_brand_but_remote_names_stay_stable() -> None:
    builder = _text("packaging/linux/build_flatpak.sh")
    remote_writer = _text("packaging/release/generate_flatpak_remote_files.py")

    assert 'REPO_TITLE="${CHEMUSON_FLATPAK_REPO_TITLE:-ChemUSON (${BRANCH})}"' in builder
    assert 'default="ChemUSON"' in remote_writer
    assert 'REPO_REMOTE_NAME="${CHEMUSON_FLATPAK_REMOTE_NAME:-chemuson-${BRANCH}}"' in builder
    assert 'REPO_CONFIG_BASENAME="${CHEMUSON_FLATPAK_CONFIG_BASENAME:-Chemuson}"' in builder
