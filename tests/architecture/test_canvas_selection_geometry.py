"""Contratos AST de la geometría extraída de la selección del canvas."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "src" / "chemuson" / "gui" / "canvas" / "canvas_selection.py"
GEOMETRY = ROOT / "src" / "chemuson" / "gui" / "canvas" / "selection_geometry.py"
CATALOG = ROOT / "architecture" / "modules.yml"

ALIASES = {
    "_signed_angle_delta_deg": "signed_angle_delta_deg",
    "_rotate_scene_point": "rotate_scene_point",
    "_optional_float_equal": "optional_float_equal",
    "_point_equal": "point_equal",
    "_scale_point_from_anchor": "scale_point_from_anchor",
    "_normalize_label_scale": "normalize_label_scale",
}
FUNCTIONS = set(ALIASES.values()) | {"normalize_custom_stroke"}


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _selection_class() -> ast.ClassDef:
    for node in _tree(SELECTION).body:
        if isinstance(node, ast.ClassDef) and node.name == "CanvasSelectionMixin":
            return node
    raise AssertionError("CanvasSelectionMixin no existe")


def test_selection_geometry_has_one_internal_owner() -> None:
    geometry_functions = {
        node.name
        for node in _tree(GEOMETRY).body
        if isinstance(node, ast.FunctionDef)
    }
    selection_methods = {
        node.name
        for node in _selection_class().body
        if isinstance(node, ast.FunctionDef)
    }

    assert FUNCTIONS <= geometry_functions
    assert set(ALIASES).isdisjoint(selection_methods)
    assert "_normalize_custom_stroke" in selection_methods


def test_selection_mixin_retains_aliases_and_stroke_wrapper() -> None:
    aliases: dict[str, str] = {}
    methods: dict[str, ast.FunctionDef] = {}
    for node in _selection_class().body:
        if isinstance(node, ast.FunctionDef):
            methods[node.name] = node
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
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

    assert aliases == ALIASES

    wrapper = methods["_normalize_custom_stroke"]
    calls = [
        node
        for node in ast.walk(wrapper)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "normalize_custom_stroke"
    ]
    assert len(calls) == 1
    assert len(calls[0].args) == 2


def test_selection_geometry_imports_only_math_and_qpointf() -> None:
    imports: set[str] = set()
    imported_names: dict[str, set[str]] = {}
    for node in _tree(GEOMETRY).body:
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            imports.add(module)
            imported_names[module] = {alias.name for alias in node.names}

    assert imports <= {"__future__", "math", "PyQt6.QtCore"}
    assert imported_names.get("PyQt6.QtCore") == {"QPointF"}


def test_selection_geometry_is_internal_to_m09() -> None:
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m09 = next(module for module in catalog["modules"] if module["id"] == "M09")

    assert "src/chemuson/gui/canvas/" in m09["paths"]
    assert "selection_geometry" in m09["internal_api"]
    assert "tests/architecture/test_canvas_selection_geometry.py" in m09["tests"]
    assert "tests/test_canvas_selection_geometry.py" in m09["tests"]
    assert not (FUNCTIONS & set(m09["public_api"]))
