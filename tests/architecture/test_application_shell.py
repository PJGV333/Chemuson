"""Contratos AST del application shell de M08, sin importar PyQt6."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MAIN_WINDOW_PATH = REPO_ROOT / "src" / "chemuson" / "gui" / "main_window.py"
SHELL_PATH = REPO_ROOT / "src" / "chemuson" / "gui" / "shell"
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _window_init() -> ast.FunctionDef:
    for node in _tree(MAIN_WINDOW_PATH).body:
        if isinstance(node, ast.ClassDef) and node.name == "ChemusonWindow":
            for member in node.body:
                if isinstance(member, ast.FunctionDef) and member.name == "__init__":
                    return member
    raise AssertionError("ChemusonWindow.__init__ no existe")


def test_main_window_delegates_shell_assembly_exactly_once() -> None:
    calls = [
        node
        for node in ast.walk(_window_init())
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "assemble_application_shell"
    ]
    assert len(calls) == 1
    assert len(calls[0].args) == 1
    assert isinstance(calls[0].args[0], ast.Name)
    assert calls[0].args[0].id == "self"


def test_constructor_no_longer_assembles_shell_regions() -> None:
    forbidden_attributes = {
        "tabs",
        "templates_dock",
        "inspector_dock",
        "validation_dock",
        "chemical_properties_dock",
        "spectroscopy_dock",
        "compchem_dock",
        "appearance_dock",
        "toolbar",
        "symbols_toolbar",
        "text_toolbar",
        "_total_charge_label",
        "_iupac_name_label",
    }
    assigned = {
        node.attr
        for node in ast.walk(_window_init())
        if isinstance(node, ast.Attribute)
        and isinstance(node.ctx, ast.Store)
        and isinstance(node.value, ast.Name)
        and node.value.id == "self"
    }
    assert assigned.isdisjoint(forbidden_attributes)


def test_shell_does_not_import_main_window() -> None:
    for path in SHELL_PATH.glob("*.py"):
        for node in ast.walk(_tree(path)):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            assert "chemuson.gui.main_window" not in names


def test_shell_is_owned_exclusively_by_m08() -> None:
    catalog = yaml.safe_load(CATALOG_PATH.read_text(encoding="utf-8"))
    owners = [
        module["id"]
        for module in catalog["modules"]
        if "src/chemuson/gui/shell/" in module["paths"]
    ]
    assert owners == ["M08"]
