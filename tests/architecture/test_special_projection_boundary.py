from __future__ import annotations

import ast
import os
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
PROJECTION_DIR = ROOT / "src/chemuson/gui/projections"


def test_projection_geometry_is_cataloged_as_internal_m08_code() -> None:
    catalog = yaml.safe_load((ROOT / "architecture/modules.yml").read_text(encoding="utf-8"))
    gui = next(module for module in catalog["modules"] if module["id"] == "M08")

    assert "src/chemuson/gui/projections/" in gui["paths"]
    assert "projections" in gui["internal_api"]
    assert "M00" in gui["current_dependencies"]


def test_projection_modules_do_not_import_qt_or_heavy_gui_renderer() -> None:
    forbidden_prefixes = (
        "PyQt6",
        "chemuson.gui.canvas",
        "chemuson.gui.dialogs",
        "chemuson.gui.items",
        "chemuson.gui.style",
    )
    for path in PROJECTION_DIR.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
            else:
                continue
            assert not any(
                module == prefix or module.startswith(prefix + ".")
                for module in modules
                for prefix in forbidden_prefixes
            ), f"heavy renderer dependency in {path}: {modules}"


def test_projection_package_import_does_not_load_qt_or_gui_widgets() -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    source = (
        "import sys; import chemuson.gui.projections; "
        "forbidden = [name for name in sys.modules if name == 'PyQt6' "
        "or name.startswith('PyQt6.') or name.startswith('chemuson.gui.canvas') "
        "or name.startswith('chemuson.gui.dialogs') or name == 'chemuson.gui.items']; "
        "assert not forbidden, forbidden"
    )
    result = subprocess.run(
        [sys.executable, "-c", source],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
