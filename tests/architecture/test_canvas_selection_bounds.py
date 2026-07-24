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
    """Prohíbe imports por prefijo: QtGui, QtWidgets, commands, dialogs,
    controllers, items (gui.), modelos de chemio, y otros."""
    imports: list[str] = []
    for node in _tree(BOUNDS).body:
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            imports.append(module)

    forbidden_prefixes = (
        "PyQt6.QtGui",
        "PyQt6.QtWidgets",
        "chemuson.gui.commands",
        "chemuson.gui.dialogs",
        "chemuson.gui.controllers",
        "chemuson.gui.items",
        "chemuson.gui.canvas",
        "chemuson.core.model",
        "chemuson.chemio",
    )
    violated = [
        imp for imp in imports
        if any(imp == fp or imp.startswith(fp + ".") or imp.startswith(fp + "/") for fp in forbidden_prefixes)
    ]
    assert violated == [], (
        f"selection_bounds importa módulos prohibidos por prefijo: {violated}"
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
    assert "extend" not in func_names, "extend no debe ser método de nivel de clase"
    assert "extend_atom_bounds" not in func_names


def test_selected_items_bbox_no_nested_extend_functions() -> None:
    """_selected_items_bbox no debe contener funciones anidadas extend/_extend
    ni extend_atom_bounds/_extend_atom_bounds. Debe delegar a selection_bounds."""
    cls = _selection_class()
    method = None
    for node in cls.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_selected_items_bbox":
            method = node
            break
    assert method is not None, "_selected_items_bbox no existe"

    # Buscar funciones anidadas dentro del cuerpo
    nested_funcs: list[str] = []
    for node in ast.walk(method):
        if isinstance(node, ast.FunctionDef) and node is not method:
            nested_funcs.append(node.name)

    forbidden_nested = [n for n in nested_funcs if n in (
        "extend", "_extend", "extend_atom_bounds", "_extend_atom_bounds"
    )]
    assert forbidden_nested == [], (
        f"_selected_items_bbox contiene funciones anidadas prohibidas: {forbidden_nested}"
    )


def test_targets_bbox_no_nested_extend_functions() -> None:
    """_targets_bbox no debe contener funciones anidadas extend/_extend
    ni extend_atom_bounds/_extend_atom_bounds. Debe delegar a selection_bounds."""
    cls = _selection_class()
    method = None
    for node in cls.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_targets_bbox":
            method = node
            break
    assert method is not None, "_targets_bbox no existe"

    nested_funcs: list[str] = []
    for node in ast.walk(method):
        if isinstance(node, ast.FunctionDef) and node is not method:
            nested_funcs.append(node.name)

    forbidden_nested = [n for n in nested_funcs if n in (
        "extend", "_extend", "extend_atom_bounds", "_extend_atom_bounds"
    )]
    assert forbidden_nested == [], (
        f"_targets_bbox contiene funciones anidadas prohibidas: {forbidden_nested}"
    )


def test_selected_atom_ids_for_transform_calls_resolve() -> None:
    """_selected_atom_ids_for_transform debe llamar a resolve_selected_atom_ids."""
    cls = _selection_class()
    method = None
    for node in cls.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_selected_atom_ids_for_transform":
            method = node
            break
    assert method is not None, "_selected_atom_ids_for_transform no existe"

    # Buscar llamadas a resolve_selected_atom_ids
    found = False
    for node in ast.walk(method):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                if node.func.attr == "resolve_selected_atom_ids":
                    found = True
            elif isinstance(node.func, ast.Name):
                if node.func.id == "resolve_selected_atom_ids":
                    found = True
    assert found, "_selected_atom_ids_for_transform no llama a resolve_selected_atom_ids"


def test_selected_items_bbox_calls_selection_bounds() -> None:
    """_selected_items_bbox debe llamar a selection_bounds."""
    cls = _selection_class()
    method = None
    for node in cls.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_selected_items_bbox":
            method = node
            break
    assert method is not None, "_selected_items_bbox no existe"

    found = False
    for node in ast.walk(method):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                if node.func.attr == "selection_bounds":
                    found = True
            elif isinstance(node.func, ast.Name):
                if node.func.id == "selection_bounds":
                    found = True
    assert found, "_selected_items_bbox no llama a selection_bounds"


def test_targets_bbox_calls_selection_bounds() -> None:
    """_targets_bbox debe llamar a selection_bounds."""
    cls = _selection_class()
    method = None
    for node in cls.body:
        if isinstance(node, ast.FunctionDef) and node.name == "_targets_bbox":
            method = node
            break
    assert method is not None, "_targets_bbox no existe"

    found = False
    for node in ast.walk(method):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                if node.func.attr == "selection_bounds":
                    found = True
            elif isinstance(node.func, ast.Name):
                if node.func.id == "selection_bounds":
                    found = True
    assert found, "_targets_bbox no llama a selection_bounds"


def test_selection_bounds_no_try_except_attribute_error() -> None:
    """selection_bounds.py no debe contener try/except AttributeError o
    try/except RuntimeError que oculten errores."""
    tree = _tree(BOUNDS)
    for node in ast.walk(tree):
        if isinstance(node, ast.Try):
            for handler in node.handlers:
                if handler.type is not None:
                    exc_names = []
                    if isinstance(handler.type, ast.Name):
                        exc_names.append(handler.type.id)
                    elif isinstance(handler.type, ast.Tuple):
                        for elt in handler.type.elts:
                            if isinstance(elt, ast.Name):
                                exc_names.append(elt.id)
                    for exc in exc_names:
                        assert exc not in (
                            "AttributeError", "RuntimeError"
                        ), f"try/except {exc} encontrado en selection_bounds.py"


def test_selection_bounds_no_any_typing() -> None:
    """selection_bounds.py no debe importar typing.Any."""
    tree = _tree(BOUNDS)
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            if node.module == "typing":
                names = {alias.name for alias in node.names}
                assert "Any" not in names, (
                    "selection_bounds.py importa typing.Any"
                )


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

    assert m09["temporary_exceptions"] == []


def test_canvas_selection_delegates_to_selection_bounds() -> None:
    """canvas_selection.py debe importar desde selection_bounds."""
    tree = _tree(SELECTION)
    import_found = False
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if "selection_bounds" in module:
                import_found = True
                names = {alias.name for alias in node.names}
                assert "resolve_selected_atom_ids" in names, "Falta resolve_selected_atom_ids"
                assert "selection_bounds" in names, "Falta selection_bounds"
                break
    assert import_found, "canvas_selection.py no importa desde selection_bounds"


def test_resolve_selected_atom_ids_uses_model_bonds_in_get_bond() -> None:
    """resolve_selected_atom_ids debe usar el patrón
    'if bond_id in model.bonds: bond = model.get_bond(bond_id)'."""
    tree = _tree(BOUNDS)
    func = None
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "resolve_selected_atom_ids":
            func = node
            break
    assert func is not None, "resolve_selected_atom_ids no existe"

    # Verificar que hay un 'if' con 'in' y luego 'get_bond'
    has_in_check = False
    has_get_bond = False
    for node in ast.walk(func):
        if isinstance(node, ast.Compare):
            for op in node.ops:
                if isinstance(op, ast.In):
                    has_in_check = True
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute) and node.func.attr == "get_bond":
                has_get_bond = True

    assert has_in_check, "resolve_selected_atom_ids no verifica 'in model.bonds'"
    assert has_get_bond, "resolve_selected_atom_ids no llama a 'model.get_bond()'"
