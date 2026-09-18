"""AST contracts for the extracted selection hit-testing policy."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "src/chemuson/gui/canvas/canvas_selection_input.py"
HIT = ROOT / "src/chemuson/gui/editor2d/selection_hit_testing.py"
CATALOG = ROOT / "architecture/modules.yml"


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _class(path: Path, name: str) -> ast.ClassDef:
    return next(node for node in _tree(path).body if isinstance(node, ast.ClassDef) and node.name == name)


def test_hit_testing_module_exists_and_owns_queries() -> None:
    assert HIT.exists()
    names = {
        node.name
        for node in _tree(HIT).body
        if isinstance(node, ast.FunctionDef)
    }
    assert names >= {
        "semantic_diagram_parent",
        "get_item_at",
        "resolve_click_item",
        "selected_annotation_item_at",
    }


def test_hit_testing_has_no_forbidden_imports_or_mutations() -> None:
    imports: list[str] = []
    for node in _tree(HIT).body:
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    forbidden = (
        "chemuson.gui.canvas.canvas_selection_input",
        "chemuson.gui.canvas.canvas_view",
        "chemuson.gui.commands",
        "chemuson.gui.controllers",
        "chemuson.gui.dialogs",
        "chemuson.core",
    )
    assert not any(module == prefix or module.startswith(prefix + ".") for module in imports for prefix in forbidden)
    mutation_names = {
        "setSelected", "setVisible", "setPos", "addItem", "removeItem",
        "clearSelection", "push", "update", "undo", "redo",
    }
    calls = {
        node.func.attr
        for node in ast.walk(_tree(HIT))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    assert calls.isdisjoint(mutation_names)


def test_input_mixin_keeps_wrappers_and_pick_consumer() -> None:
    cls = _class(INPUT, "CanvasSelectionInputMixin")
    methods = {
        node.name: node
        for node in cls.body
        if isinstance(node, ast.FunctionDef)
    }
    for name in (
        "_semantic_diagram_parent",
        "_get_item_at",
        "_resolve_click_item",
        "_selected_annotation_item_at",
    ):
        assert name in methods
    assert any(
        isinstance(node, ast.Attribute)
        and node.attr == "_pick_hover_target"
        for node in ast.walk(methods["_resolve_click_item"])
    )
    assert any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id in {
            "semantic_diagram_parent",
            "get_item_at",
            "resolve_click_item",
            "selected_annotation_item_at",
        }
        for method in methods.values()
        for node in ast.walk(method)
    )


def test_mouse_press_event_remains_in_input_mixin() -> None:
    cls = _class(INPUT, "CanvasSelectionInputMixin")
    assert any(
        isinstance(node, ast.FunctionDef) and node.name == "mousePressEvent"
        for node in cls.body
    )


def test_catalog_records_hit_testing_under_m20_without_debt() -> None:
    import yaml

    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m20 = next(module for module in catalog["modules"] if module["id"] == "M20")
    assert "selection_hit_testing" in m20["internal_api"]
    assert "tests/test_canvas_selection_hit_testing.py" in m20["tests"]
    assert "tests/architecture/test_canvas_selection_hit_testing.py" in m20["tests"]
    assert m20["temporary_exceptions"] == []
    assert m20["circular_dependencies"] == []
