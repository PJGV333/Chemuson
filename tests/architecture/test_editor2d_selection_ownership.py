"""Architecture contracts for the M20 editor2d selection helper package."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "architecture" / "modules.yml"
HELPERS = (
    "selection_geometry",
    "selection_bounds",
    "selection_hit_testing",
    "selection_overlay",
    "selection_clipboard",
)


def _catalog() -> list[dict]:
    return yaml.safe_load(CATALOG.read_text(encoding="utf-8"))["modules"]


def _module(module_id: str) -> dict:
    return next(module for module in _catalog() if module["id"] == module_id)


def test_m20_canonical_helpers_exist() -> None:
    for helper in HELPERS:
        path = ROOT / "src" / "chemuson" / "gui" / "editor2d" / f"{helper}.py"
        assert path.exists(), f"Missing canonical helper: {path}"


def test_legacy_modules_are_import_only_shims() -> None:
    for helper in HELPERS:
        path = ROOT / "src" / "chemuson" / "gui" / "canvas" / f"{helper}.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        definitions = [
            node
            for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        ]
        assert not definitions, f"Legacy helper owns definitions: {path}"
        imports = [node for node in tree.body if isinstance(node, ast.ImportFrom)]
        assert any(
            node.module == f"chemuson.gui.editor2d.{helper}" for node in imports
        ), f"Legacy helper does not re-export M20: {path}"


def test_m20_catalog_owns_helpers_and_tests() -> None:
    m20 = _module("M20")
    m09 = _module("M09")
    assert m20["paths"] == ["src/chemuson/gui/editor2d/"]
    assert set(m20["internal_api"]) == set(HELPERS)
    assert "src/chemuson/gui/editor2d/" not in m09["paths"]
    assert all(helper not in m09["internal_api"] for helper in HELPERS)
    assert "M20" in m09["current_dependencies"]
    assert "M20" in m09["target_dependencies"]
    assert m20["current_dependencies"] == []
    assert m20["target_dependencies"] == []
    assert m20["temporary_exceptions"] == []
    assert m20["circular_dependencies"] == []
    for helper in HELPERS:
        assert any(helper in test for test in m20["tests"])


def test_canvas_selection_imports_canonical_editor2d_helpers() -> None:
    expected = {
        "canvas_selection.py": {
            "chemuson.gui.editor2d.selection_geometry",
            "chemuson.gui.editor2d.selection_bounds",
            "chemuson.gui.editor2d.selection_overlay",
            "chemuson.gui.editor2d.selection_clipboard",
        },
        "canvas_selection_input.py": {
            "chemuson.gui.editor2d.selection_hit_testing",
        },
    }
    for filename, modules in expected.items():
        tree = ast.parse(
            (ROOT / "src" / "chemuson" / "gui" / "canvas" / filename).read_text(
                encoding="utf-8"
            )
        )
        imported = {
            node.module
            for node in tree.body
            if isinstance(node, ast.ImportFrom) and node.module
        }
        assert modules <= imported


def test_all_five_helpers_have_unique_canonical_owner() -> None:
    catalog_text = CATALOG.read_text(encoding="utf-8")
    assert catalog_text.count("src/chemuson/gui/editor2d/") == 1
    m20 = _module("M20")
    assert set(m20["internal_api"]) == set(HELPERS)
