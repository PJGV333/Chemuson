from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def _env(config_home: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["XDG_CONFIG_HOME"] = str(config_home)
    paths = [str(ROOT), str(ROOT / "src")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def test_lab_window_starts_and_closes_in_bounded_offscreen_process(tmp_path) -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "tools/visual/depiction_lab.py"), "--smoke"],
        cwd=ROOT,
        env=_env(tmp_path / "smoke-config"),
        capture_output=True,
        text=True,
        timeout=30.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_ui_initial_centering_and_both_viewports_stay_visible_after_resize_and_transform(tmp_path) -> None:
    probe = ROOT / "tests/probes/depiction_lab_viewport.py"
    result = subprocess.run(
        [sys.executable, str(probe)],
        cwd=ROOT,
        env=_env(tmp_path / "viewport-config"),
        capture_output=True,
        text=True,
        timeout=90.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '"observations": 12' in result.stdout
    assert '"closed": true' in result.stdout


def test_ui_parameter_updates_case_selection_zoom_theme_and_profile_lifecycle(tmp_path) -> None:
    probe = ROOT / "tests/probes/depiction_lab_ui.py"
    result = subprocess.run(
        [sys.executable, str(probe)],
        cwd=ROOT,
        env=_env(tmp_path / "ui-config"),
        capture_output=True,
        text=True,
        timeout=75.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '"cases": 12' in result.stdout
    assert '"closed": true' in result.stdout
