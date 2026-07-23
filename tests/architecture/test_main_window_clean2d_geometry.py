"""Contratos AST de la geometría Clean2D extraída de la ventana."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
MAIN_WINDOW = ROOT / "src" / "chemuson" / "gui" / "main_window.py"
GEOMETRY = ROOT / "src" / "chemuson" / "gui" / "clean2d_geometry.py"
CATALOG = ROOT / "architecture" / "modules.yml"

HELPERS = {
    "_coords_center": "coords_center",
    "_average_bond_length": "average_bond_length",
    "_rescale_coords_to_bond_length": "rescale_coords_to_bond_length",
    "_align_coords_to_reference": "align_coords_to_reference",
    "_project_missing_hydrogen_coords": "project_missing_hydrogen_coords",
}


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _window_class() -> ast.ClassDef:
    for node in _tree(MAIN_WINDOW).body:
        if isinstance(node, ast.ClassDef) and node.name == "ChemusonWindow":
            return node
    raise AssertionError("ChemusonWindow no existe")


def test_geometry_helpers_have_one_pure_owner() -> None:
    geometry_functions = {
        node.name
        for node in _tree(GEOMETRY).body
        if isinstance(node, ast.FunctionDef)
    }
    window_methods = {
        node.name
        for node in _window_class().body
        if isinstance(node, ast.FunctionDef)
    }

    assert set(HELPERS.values()) <= geometry_functions
    assert set(HELPERS).isdisjoint(window_methods)


def test_main_window_retains_static_private_aliases() -> None:
    aliases: dict[str, str] = {}
    for node in _window_class().body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        value = node.value
        if (
            isinstance(target, ast.Name)
            and isinstance(value, ast.Call)
            and isinstance(value.func, ast.Name)
            and value.func.id == "staticmethod"
            and len(value.args) == 1
            and isinstance(value.args[0], ast.Name)
        ):
            aliases[target.id] = value.args[0].id

    assert aliases == HELPERS


def test_geometry_module_imports_only_standard_library() -> None:
    imports: set[str] = set()
    for node in _tree(GEOMETRY).body:
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module or "")

    assert imports <= {"__future__", "math"}


def test_clean2d_geometry_is_internal_to_m08() -> None:
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m08 = next(module for module in catalog["modules"] if module["id"] == "M08")

    assert "src/chemuson/gui/clean2d_geometry.py" in m08["paths"]
    assert "clean2d_geometry" in m08["internal_api"]
    assert "tests/architecture/test_main_window_clean2d_geometry.py" in m08["tests"]
    assert m08["public_api"] == ["ChemusonWindow"]
