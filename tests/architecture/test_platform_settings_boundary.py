from __future__ import annotations

import ast
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
PLATFORM_PATH = ROOT / "src" / "chemuson" / "platform"


def _catalog() -> dict:
    return yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text())


def _module(module_id: str) -> dict:
    return next(item for item in _catalog()["modules"] if item["id"] == module_id)


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    result: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            result.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            result.add(node.module or "")
    return result


def test_platform_settings_has_no_gui_or_widget_dependencies() -> None:
    for path in PLATFORM_PATH.glob("*.py"):
        imports = _imports(path)
        assert not any(name.startswith("chemuson.gui") for name in imports)
        assert "PyQt6.QtWidgets" not in imports


def test_platform_resources_is_canonical_and_legacy_helper_is_a_shim() -> None:
    canonical = (PLATFORM_PATH / "resources.py").read_text()
    legacy = (ROOT / "src" / "chemuson" / "utils" / "resources.py").read_text()
    assert "def open_resource_path" in canonical
    assert "from chemuson.platform.resources import open_resource_path" in legacy
    assert "def open_resource_path" not in legacy


def test_platform_catalog_lists_explicit_files() -> None:
    platform = _module("M21")
    assert set(platform["paths"]) == {
        "src/chemuson/platform/__init__.py",
        "src/chemuson/platform/settings.py",
        "src/chemuson/platform/resources.py",
    }


def test_platform_catalog_and_dependency_direction_are_explicit() -> None:
    platform = _module("M21")
    gui = _module("M08")
    utils = _module("M15")
    assert platform["name"] == "platform.settings"
    assert platform["paths"] == [
        "src/chemuson/platform/__init__.py",
        "src/chemuson/platform/settings.py",
        "src/chemuson/platform/resources.py",
    ]
    assert platform["current_dependencies"] == []
    assert platform["target_dependencies"] == []
    assert platform["temporary_exceptions"] == []
    assert platform["circular_dependencies"] == []
    assert "M08" in platform["forbidden_dependencies"]
    assert "M21" in gui["current_dependencies"]
    assert "M21" in gui["target_dependencies"]
    assert "M21" in utils["current_dependencies"]
    assert "M21" in utils["target_dependencies"]
