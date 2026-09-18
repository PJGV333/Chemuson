"""AST contracts for extracted selection overlay calculations."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "src" / "chemuson" / "gui" / "canvas" / "canvas_selection.py"
OVERLAY = ROOT / "src" / "chemuson" / "gui" / "editor2d" / "selection" / "selection_overlay.py"
CATALOG = ROOT / "architecture" / "modules.yml"


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _selection_class() -> ast.ClassDef:
    return next(
        node for node in _tree(SELECTION).body
        if isinstance(node, ast.ClassDef) and node.name == "CanvasSelectionMixin"
    )


def test_overlay_module_owns_calculation_functions() -> None:
    names = {
        node.name for node in _tree(OVERLAY).body if isinstance(node, ast.FunctionDef)
    }
    assert names >= {
        "padded_selection_bbox",
        "offset_scene_point",
        "selection_handle_scene_positions",
        "handle_item_distance_sq",
        "handle_item_hit_radius",
        "selection_handle_hit_kind",
    }


def test_overlay_module_has_no_forbidden_imports_or_scene_mutations() -> None:
    imports: list[str] = []
    for node in _tree(OVERLAY).body:
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    forbidden_imports = (
        "PyQt6.QtGui", "PyQt6.QtWidgets", "chemuson.gui.canvas",
        "chemuson.gui.commands", "chemuson.gui.controllers", "chemuson.gui.dialogs",
    )
    assert not any(module == prefix or module.startswith(prefix + ".") for module in imports for prefix in forbidden_imports)
    mutation_names = {
        "setSelected", "setVisible", "setPos", "setRect", "addItem", "removeItem",
        "clearSelection", "push", "update", "undo", "redo",
    }
    calls = {
        node.func.attr
        for node in ast.walk(_tree(OVERLAY))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    assert calls.isdisjoint(mutation_names)


def test_canvas_keeps_scene_mutation_and_private_wrappers() -> None:
    names = {
        node.name for node in _selection_class().body if isinstance(node, ast.FunctionDef)
    }
    for name in (
        "_ensure_selection_overlay", "_apply_selection_overlay_bbox",
        "_update_selection_overlay", "_update_drag_selection_overlay",
        "_hit_selection_handle", "_hit_selection_move_handle",
        "_hit_selection_scale_handle", "_hit_handle_item",
        "_handle_item_distance_sq", "_handle_item_hit_radius",
        "_selection_handle_hit_kind",
    ):
        assert name in names


def test_canvas_delegates_overlay_calculations() -> None:
    names = {
        "padded_selection_bbox", "offset_scene_point", "selection_handle_scene_positions",
        "handle_item_distance_sq", "handle_item_hit_radius", "selection_handle_hit_kind",
    }
    calls = {
        node.func.id
        for node in ast.walk(_tree(SELECTION))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert names <= calls


def test_catalog_records_overlay_module_and_tests() -> None:
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m20 = next(module for module in catalog["modules"] if module["id"] == "M20")
    assert "selection_overlay" in m20["internal_api"]
    assert "tests/test_canvas_selection_overlay.py" in m20["tests"]
    assert "tests/architecture/test_canvas_selection_overlay.py" in m20["tests"]
    assert m20["temporary_exceptions"] == []
    assert m20["circular_dependencies"] == []
