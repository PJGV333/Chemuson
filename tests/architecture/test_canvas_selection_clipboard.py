"""AST contracts for the extracted clipboard policy module."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "src" / "chemuson" / "gui" / "canvas" / "canvas_selection.py"
CLIPBOARD = ROOT / "src" / "chemuson" / "gui" / "editor2d" / "selection_clipboard.py"
CATALOG = ROOT / "architecture" / "modules.yml"


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _selection_class() -> ast.ClassDef:
    return next(node for node in _tree(SELECTION).body if isinstance(node, ast.ClassDef) and node.name == "CanvasSelectionMixin")


def test_clipboard_module_owns_policy_and_codec_functions() -> None:
    names = {node.name for node in _tree(CLIPBOARD).body if isinstance(node, ast.FunctionDef)}
    assert names >= {
        "mime_has_pasteable_format", "encode_selection_payload", "decode_selection_payload",
        "is_large_clipboard_structure", "bond_copy_priority", "unique_bonds_for_copy",
    }


def test_clipboard_module_has_no_commands_undo_or_scene_mutations() -> None:
    imports: list[str] = []
    for node in _tree(CLIPBOARD).body:
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    forbidden = (
        "chemuson.gui.commands", "chemuson.gui.controllers", "chemuson.gui.dialogs",
        "chemuson.gui.canvas.canvas_selection", "PyQt6.QtWidgets",
    )
    assert not any(module == prefix or module.startswith(prefix + ".") for module in imports for prefix in forbidden)
    mutation_names = {"setSelected", "setVisible", "setPos", "addItem", "removeItem", "clearSelection", "push", "undo", "redo", "update"}
    calls = {node.func.attr for node in ast.walk(_tree(CLIPBOARD)) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)}
    assert calls.isdisjoint(mutation_names)


def test_clipboard_wrappers_and_paste_coordinator_remain() -> None:
    names = {node.name for node in _selection_class().body if isinstance(node, ast.FunctionDef)}
    for name in ("has_copyable_selection", "can_paste_from_clipboard", "copy_to_clipboard", "paste_from_clipboard", "_paste_selection_payload"):
        assert name in names
    assert "_build_selection_graph" in names


def test_clipboard_codec_calls_and_mime_wrappers_are_present() -> None:
    calls = {node.func.id for node in ast.walk(_tree(SELECTION)) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}
    assert {"encode_selection_payload", "decode_selection_payload", "mime_has_pasteable_format"} <= calls


def test_catalog_records_clipboard_policy_module_and_tests() -> None:
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m20 = next(module for module in catalog["modules"] if module["id"] == "M20")
    assert "selection_clipboard" in m20["internal_api"]
    assert "tests/test_canvas_selection_clipboard_policy.py" in m20["tests"]
    assert "tests/architecture/test_canvas_selection_clipboard.py" in m20["tests"]
    assert m20["temporary_exceptions"] == []
    assert m20["circular_dependencies"] == []
