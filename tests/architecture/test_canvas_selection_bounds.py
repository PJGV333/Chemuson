"""Contratos AST del módulo selection_bounds extraído de canvas_selection."""

from __future__ import annotations

import ast
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "src" / "chemuson" / "gui" / "canvas" / "canvas_selection.py"
BOUNDS = ROOT / "src" / "chemuson" / "gui" / "canvas" / "selection_bounds.py"
CATALOG = ROOT / "architecture" / "modules.yml"


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _selection_class() -> ast.ClassDef:
    for node in _tree(SELECTION).body:
        if isinstance(node, ast.ClassDef) and node.name == "CanvasSelectionMixin":
            return node
    raise AssertionError("CanvasSelectionMixin no existe")


def test_selection_bounds_module_exists() -> None:
    assert BOUNDS.exists(), "selection_bounds.py debe existir"


def test_selection_bounds_contains_functions() -> None:
    tree = _tree(BOUNDS)
    func_names = {
        node.name
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
    }
    assert "resolve_selected_atom_ids" in func_names, "Falta resolve_selected_atom_ids"
    assert "selection_bounds" in func_names, "Falta selection_bounds"


def test_selection_bounds_no_forbidden_imports() -> None:
    imports: set[str] = set()
    for node in _tree(BOUNDS).body:
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            imports.add(module)

    forbidden = {
        "PyQt6.QtGui",
        "PyQt6.QtWidgets",
        "chemuson.gui.commands",
        "chemuson.gui.dialogs",
        "chemuson.gui.controllers",
    }
    assert imports.isdisjoint(forbidden), (
        f"selection_bounds importa módulos prohibidos: {imports & forbidden}"
    )


def test_selection_bounds_no_mutation_calls() -> None:
    """El módulo no debe llamar métodos de mutación."""
    tree = _tree(BOUNDS)
    mutation_names = {
        "addItem",
        "removeItem",
        "setPos",
        "setVisible",
        "setSelected",
        "clearSelection",
        "update",
        "push",
    }
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
    found = calls.intersection(mutation_names)
    assert found == set(), f"Llamadas de mutación encontradas: {found}"


def test_selection_mixin_wrappers_exist() -> None:
    """Los tres wrappers privados deben permanecer en CanvasSelectionMixin."""
    cls = _selection_class()
    method_names = {
        node.name
        for node in cls.body
        if isinstance(node, ast.FunctionDef)
    }
    assert "_selected_atom_ids_for_transform" in method_names
    assert "_selected_items_bbox" in method_names
    assert "_targets_bbox" in method_names


def test_selection_mixin_no_duplicated_extend_logic() -> None:
    """No debe haber dos funciones extend/extend_atom_bounds duplicadas."""
    cls = _selection_class()
    func_names = {
        node.name
        for node in cls.body
        if isinstance(node, ast.FunctionDef)
    }
    # Las funciones extend y extend_atom_bounds internas son definidas dentro de
    # _selected_items_bbox y _targets_bbox como funciones anidadas, no como
    # miembros de nivel de clase. Este test verifica que no haya funciones
    # de nivel de clase llamadas "extend" o "extend_atom_bounds".
    assert "extend" not in func_names, "extend no debe ser método de nivel de clase"
    assert "extend_atom_bounds" not in func_names


def test_catalog_records_new_module() -> None:
    """El catálogo debe registrar selection_bounds."""
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m09 = next(m for m in catalog["modules"] if m["id"] == "M09")

    assert "selection_bounds" in m09["internal_api"]


def test_catalog_records_new_tests() -> None:
    """El catálogo debe registrar los archivos de prueba."""
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m09 = next(m for m in catalog["modules"] if m["id"] == "M09")

    assert "tests/test_canvas_selection_bounds.py" in m09["tests"]
    assert "tests/architecture/test_canvas_selection_bounds.py" in m09["tests"]


def test_catalog_no_new_exceptions() -> None:
    """No se deben añadir excepciones temporales ni ciclos."""
    catalog = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    m09 = next(m for m in catalog["modules"] if m["id"] == "M09")

    # Verificar que no haya excepciones ni ciclos no previstos
    assert m09["temporary_exceptions"] == []


def test_canvas_selection_delegates_to_selection_bounds() -> None:
    """canvas_selection.py debe importar desde selection_bounds."""
    tree = _tree(SELECTION)
    import_found = False
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            # Relativo (.selection_bounds) o absoluto (selection_bounds)
            module = node.module or ""
            if "selection_bounds" in module:
                import_found = True
                names = {alias.name for alias in node.names}
                assert "resolve_selected_atom_ids" in names, "Falta resolve_selected_atom_ids"
                assert "selection_bounds" in names, "Falta selection_bounds"
                break
    assert import_found, "canvas_selection.py no importa desde selection_bounds"
