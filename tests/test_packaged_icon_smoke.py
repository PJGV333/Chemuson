"""Contracts for frozen Windows/Linux SVG resources and their release gates."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import yaml

from chemuson.gui.theme.icon_smoke import run_icon_smoke_test

ROOT = Path(__file__).resolve().parent.parent


def test_source_smoke_renders_all_static_and_essential_icons() -> None:
    report = run_icon_smoke_test()

    assert report["status"] == "success"
    assert report["frozen"] == "false"
    assert report["static_svg_count"] == 69
    assert set(report["themes"]) == {"light", "dark"}
    assert report["hidpi_dpr"] == 2.0
    for theme in report["themes"].values():
        assert theme["static_icons_with_visible_pixels"] == 69
        assert theme["static_svg_count"] == 69
        assert all(
            value > 0
            for name, value in theme.items()
            if name.startswith(("essential:", "toolbar:"))
        )


def test_frozen_smoke_cli_is_gated_by_explicit_test_environment() -> None:
    environment = os.environ.copy()
    environment.pop("CHEMUSON_ICON_SMOKE_TEST", None)
    result = subprocess.run(
        [sys.executable, "-m", "chemuson", "--icon-smoke-test"],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 2
    assert "reserved for the packaging smoke workflow" in result.stderr


def test_pyinstaller_spec_explicitly_places_all_static_icons_under_package_path() -> None:
    spec = (ROOT / "chemuson.spec").read_text(encoding="utf-8")
    assert 'STATIC_ICON_DIR = PROJECT_ROOT / "src" / "chemuson" / "gui" / "theme" / "icons"' in spec
    assert 'len(STATIC_ICON_FILES) != 69' in spec
    assert '(str(path), "chemuson/gui/theme/icons") for path in STATIC_ICON_FILES' in spec
    assert 'datas = datas_c + datas_icons + datas_qt' in spec


def test_preview_and_release_gate_frozen_icon_smokes_on_both_platforms() -> None:
    preview_text = (ROOT / ".github/workflows/build-preview.yml").read_text(encoding="utf-8")
    release_text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    preview = yaml.load(preview_text, Loader=yaml.BaseLoader)
    release = yaml.load(release_text, Loader=yaml.BaseLoader)

    assert "validate_packaged_icons.py" in preview_text
    assert '--executable "dist/Chemuson.exe"' in preview_text
    assert "--executable dist/Chemuson" in preview_text
    assert "validate_packaged_icons.py" in release_text
    assert '--executable "dist/Chemuson.exe"' in release_text
    assert "--executable dist/Chemuson" in release_text
    assert "validate_packaged_icons.py" in (ROOT / "packaging/release/validate_appimage.py").read_text()

    windows_preview_steps = preview["jobs"]["build_windows"]["steps"]
    linux_preview_steps = preview["jobs"]["build_linux_appimage"]["steps"]
    assert any("validate_packaged_icons.py" in step.get("run", "") for step in windows_preview_steps)
    assert any("validate_packaged_icons.py" in step.get("run", "") for step in linux_preview_steps)

    windows_release_steps = release["jobs"]["build_windows"]["steps"]
    linux_release_steps = release["jobs"]["build_linux"]["steps"]
    assert any("validate_packaged_icons.py" in step.get("run", "") for step in windows_release_steps)
    assert any("validate_packaged_icons.py" in step.get("run", "") for step in linux_release_steps)


def test_frozen_icon_validator_requires_every_theme_and_nonblank_essential() -> None:
    validator = (ROOT / "packaging/release/validate_packaged_icons.py").read_text(encoding="utf-8")
    smoke = (ROOT / "src/chemuson/gui/theme/icon_smoke.py").read_text(encoding="utf-8")

    assert "EXPECTED_SVG_COUNT = 69" in validator
    assert '"static_icons_with_visible_pixels"' in smoke
    assert '"light", "dark"' in smoke
    assert "_tool_button_image" in smoke
    assert 'f"toolbar:{name}"' in validator
    assert "devicePixelRatio() != 2.0" in smoke
    for icon_name in (
        "pointer/select",
        "single bond",
        "aromatic ring",
        "search",
        "undo",
        "redo",
        "new document",
        "clean",
    ):
        assert icon_name in smoke
