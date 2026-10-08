"""Regression tests for Type 2 AppImage generation and validation contracts."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
RELEASE_DIR = ROOT / "packaging" / "release"
sys.path.insert(0, str(RELEASE_DIR))

import validate_appimage  # noqa: E402


def _elf_header(*, marker: bytes = b"AI\x02", machine: int = 62, elf_class: int = 2) -> bytes:
    header = bytearray(64)
    header[:8] = b"\x7fELF" + bytes((elf_class, 1, 1, 0))
    header[8:11] = marker
    header[16:18] = (2).to_bytes(2, "little")
    header[18:20] = machine.to_bytes(2, "little")
    return bytes(header)


def test_type2_header_requires_elf_x86_64_and_exact_ai02_signature(tmp_path: Path) -> None:
    appimage = tmp_path / "Chemuson.AppImage"
    payload = _elf_header() + b"squashfs payload"
    appimage.write_bytes(payload)

    assert validate_appimage.validate_type2_header(appimage) == hashlib.sha256(payload).hexdigest()

    for name, invalid in (
        ("no-elf.AppImage", b"not an ELF" + b"\\x00" * 60),
        ("wrong-magic.AppImage", _elf_header(marker=b"AI\x01") + b"payload"),
        ("wrong-arch.AppImage", _elf_header(machine=3) + b"payload"),
        ("wrong-class.AppImage", _elf_header(elf_class=1) + b"payload"),
    ):
        candidate = tmp_path / name
        candidate.write_bytes(invalid)
        with pytest.raises(ValueError):
            validate_appimage.validate_type2_header(candidate)


def test_type2_validator_rejects_extension_only_and_truncated_packages(tmp_path: Path) -> None:
    renamed = tmp_path / "renamed.AppImage"
    renamed.write_bytes(b"PyInstaller executable renamed to .AppImage")
    with pytest.raises(ValueError, match="ELF"):
        validate_appimage.validate_type2_header(renamed)

    truncated = tmp_path / "truncated.AppImage"
    truncated.write_bytes(b"\x7fELF")
    with pytest.raises(ValueError, match="ELF"):
        validate_appimage.validate_type2_header(truncated)

    wrong_suffix = tmp_path / "Chemuson.bin"
    wrong_suffix.write_bytes(_elf_header() + b"payload")
    with pytest.raises(ValueError, match="suffix"):
        validate_appimage.validate_type2_header(wrong_suffix)


def test_appimage_builder_uses_appdir_and_pinned_tool_not_copy_or_rename() -> None:
    builder = (ROOT / "packaging/linux/build_appimage.sh").read_text(encoding="utf-8")
    tool_fetcher = (ROOT / "packaging/linux/fetch_appimagetool.sh").read_text(encoding="utf-8")
    assert "usr/bin/Chemuson" in builder
    assert "AppRun" in builder
    assert '${APP_ID}.appdata.xml' in builder
    desktop_template = (ROOT / "packaging/linux/appimage/io.github.PJGV333.Chemuson.desktop.in").read_text()
    assert 'Categories=Science;Chemistry;' in desktop_template
    assert 'ARCH=x86_64 "${APPIMAGETOOL}" "${APPDIR}" "${APPIMAGE_PATH}"' in builder
    assert "--updateinformation \"${UPDATE_INFO}\"" in builder
    assert 'cp "$src" "$APPIMAGE_PATH"' not in builder
    assert 'cp "${SOURCE_BINARY}" "${APPIMAGE_PATH}"' not in builder
    assert "AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage" in tool_fetcher
    assert "98605504" in tool_fetcher
    assert "5735cc5" in tool_fetcher
    assert "b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2" in tool_fetcher
    assert "sha256sum --check" in tool_fetcher
    assert "--appimage-extract" in tool_fetcher


def test_preview_and_release_workflows_build_and_validate_same_type2_package() -> None:
    preview_text = (ROOT / ".github/workflows/build-preview.yml").read_text(encoding="utf-8")
    release_text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    preview = yaml.load(preview_text, Loader=yaml.BaseLoader)
    release = yaml.load(release_text, Loader=yaml.BaseLoader)

    assert preview["jobs"]["build_linux_appimage"]["runs-on"] == "ubuntu-22.04"
    assert release["jobs"]["build_linux"]["runs-on"] == "ubuntu-22.04"
    for workflow_text in (preview_text, release_text):
        assert "packaging/linux/build_appimage.sh" in workflow_text
        assert "packaging/release/validate_appimage.py" in workflow_text
        assert "--appimage-extract" in (ROOT / "packaging/release/validate_appimage.py").read_text()
    assert "--build-type preview" in preview_text
    assert "--build-type release" in release_text
    assert "appstream" in preview_text
    assert "appstream" in release_text
    assert "desktop-file-utils" in preview_text
    assert "desktop-file-utils" in release_text
    assert "zsyncmake" in (ROOT / "packaging/linux/build_appimage.sh").read_text()


def test_release_validator_requires_type2_and_zsync_without_changing_asset_name() -> None:
    validator = (ROOT / "packaging/release/validate_release_artifacts.py").read_text(encoding="utf-8")
    builder = (ROOT / "packaging/linux/build_appimage.sh").read_text(encoding="utf-8")
    assert "validate_type2_header(appimage_path)" in validator
    assert 'AppImage.zsync' in validator
    assert 'Chemuson-v${VERSION}-linux-x86_64.AppImage' in builder
    assert 'gh-releases-zsync|${OWNER}|${REPO}|${UPDATE_TRACK}|' in builder
    assert 'appimage_update_information' in builder


def test_appdir_requires_launcher_desktop_icon_and_appstream(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="AppDir is missing"):
        validate_appimage._validate_appdir(tmp_path, version="0.3.0-beta.1")


def test_release_gate_checks_appimage_packaging_contracts() -> None:
    release_text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    validator = (ROOT / "packaging/release/validate_appimage.py").read_text()
    assert "tests/test_appimage_validation.py" in release_text
    assert "tests/test_packaged_icon_smoke.py" in release_text
    assert "packaging/release/validate_appimage.py" in release_text
    assert "validate_packaged_icons.py" in validator
    assert "PyInstaller CArchive must contain all 69 static SVG icons" in validator
