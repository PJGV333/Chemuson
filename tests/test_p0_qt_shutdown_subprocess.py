from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _isolated_environment(config_dir: Path) -> dict[str, str]:
    environment = os.environ.copy()
    paths = [str(ROOT / "src")]
    rdkit_packages = ROOT / "chem/lib/python3.14/site-packages"
    if rdkit_packages.is_dir():
        paths.insert(0, str(rdkit_packages))
    if environment.get("PYTHONPATH"):
        paths.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(paths)
    environment["QT_QPA_PLATFORM"] = "offscreen"
    environment["XDG_CONFIG_HOME"] = str(config_dir)
    return environment


def test_file_workflow_close_and_autosave_recovery_in_bounded_processes(tmp_path) -> None:
    probe = ROOT / "tests/probes/p0_shutdown_workflow.py"
    scenarios = ("normal", "exported", "modified", "canvas-delete", "async-cleanup", "recovered")
    for scenario in scenarios:
        result = subprocess.run(
            [sys.executable, str(probe), scenario],
            capture_output=True,
            text=True,
            timeout=25.0,
            env=_isolated_environment(tmp_path / f"config-{scenario}"),
            check=False,
        )
        assert result.returncode == 0, f"{scenario}:\n{result.stdout}{result.stderr}"
        assert '"closed": true' in result.stdout, result.stdout
        assert '"window_deleted": true' in result.stdout, result.stdout


def test_existing_async_worker_shutdown_contract_in_isolated_pytest_process(tmp_path) -> None:
    pytest_executable = shutil.which("pytest")
    assert pytest_executable is not None
    result = subprocess.run(
        [pytest_executable, "-q", str(ROOT / "tests/test_gui_async_worker_shutdown.py")],
        capture_output=True,
        text=True,
        timeout=90.0,
        cwd=ROOT,
        env=_isolated_environment(tmp_path / "async-contract"),
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "2 passed" in result.stdout
