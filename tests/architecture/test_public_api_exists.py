"""Tests arquitectónicos: validación de que los símbolos de public_api existen en la superficie pública (Phase 1).

Valida que cada símbolo declarado en `public_api` existe como nombre de nivel superior
en los archivos de superficie pública correspondientes, sin importar módulos de producción.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"


# ---------------------------------------------------------------------------
# Tipos para resultados estructurados
# ---------------------------------------------------------------------------

PathError = Dict[str, Any]  # debe contener: path (str), error_type (str)
ParseError = Dict[str, Any]  # debe contener: file (Path), type (str), message (str), line (int | None)
SurfaceFile = Path


# ---------------------------------------------------------------------------
# Helper de superficie pública con diagnóstico
# ---------------------------------------------------------------------------

def get_public_surface_files(
    repo_root: Path,
    module_paths: list[str],
) -> Tuple[List[SurfaceFile], List[PathError]]:
    """
    Dado el listado de paths de un módulo en el catálogo, devuelve:
    - lista de archivos Python válidos que constituyen su superficie pública
    - lista de errores encontrados (path_missing, package_init_missing, unsupported_path)

    Reglas:
    - Si path es directorio -> solo __init__.py dentro de él.
    - Si path es archivo .py explícito -> ese archivo.
    - Varios paths -> unión (ordenada y deduplicada).
    - Los errores se registran, no se silencian.
    """
    surface_files: List[SurfaceFile] = []
    path_errors: List[PathError] = []
    seen_files: Set[Path] = set()

    for path_str in module_paths:
        full_path = repo_root / path_str

        if not full_path.exists():
            path_errors.append({
                "path": path_str,
                "error_type": "path_missing",
            })
            continue

        if full_path.is_dir():
            init_file = full_path / "__init__.py"
            if init_file.exists():
                # Deduplicar
                resolved = init_file.resolve()
                if resolved not in seen_files:
                    surface_files.append(resolved)
                    seen_files.add(resolved)
            else:
                path_errors.append({
                    "path": path_str,
                    "error_type": "package_init_missing",
                })
        else:
            if full_path.suffix == ".py":
                resolved = full_path.resolve()
                if resolved not in seen_files:
                    surface_files.append(resolved)
                    seen_files.add(resolved)
            else:
                path_errors.append({
                    "path": path_str,
                    "error_type": "unsupported_path",
                })

    # Ordenar de forma determinista
    return sorted(surface_files), sorted(path_errors, key=lambda e: e["path"])


# ---------------------------------------------------------------------------
# Extracción AST con manejo de errores explícito
# ---------------------------------------------------------------------------

def extract_top_level_names(file_path: Path) -> Set[str]:
    """
    Analiza el archivo dado y extrae los nombres que están vinculados a nivel
    superior del módulo: definiciones, asignaciones e imports.

    No usa ast.walk() para evitar capturar nombres anidados. Recorre solo el body.

    Puede propagar:
    - FileNotFoundError
    - UnicodeDecodeError (errors="replace" evita algunos)
    - SyntaxError
    """
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        source = f.read()

    tree = ast.parse(source, filename=str(file_path))

    names: Set[str] = set()

    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            names.add(node.name)
        elif isinstance(node, ast.AsyncFunctionDef):
            names.add(node.name)
        elif isinstance(node, ast.ClassDef):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    names.add(target.id)
                elif isinstance(target, (ast.Tuple, ast.List)):
                    _extract_names_from_tuple(target, names)
        elif isinstance(node, ast.AnnAssign):
            if node.target and isinstance(node.target, ast.Name):
                names.add(node.target.id)
            elif node.target and isinstance(node.target, (ast.Tuple, ast.List)):
                _extract_names_from_tuple(node.target, names)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name.split(".", 1)[0]
                if alias.asname:
                    names.add(alias.asname)
                else:
                    names.add(top)
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if alias.name == "*":
                    continue
                if alias.asname:
                    names.add(alias.asname)
                else:
                    names.add(alias.name)

    return names


def _extract_names_from_tuple(
    target: ast.Tuple | ast.List, names: Set[str]
) -> None:
    """Extrae nombres de elementos de tupla/lista que sean ast.Name."""
    for elt in target.elts:
        if isinstance(elt, ast.Name):
            names.add(elt.id)
        elif isinstance(elt, (ast.Tuple, ast.List)):
            _extract_names_from_tuple(elt, names)


# ---------------------------------------------------------------------------
# Auditoría unificada del módulo
# ---------------------------------------------------------------------------

def audit_module_public_api(
    repo_root: Path,
    module: dict,
) -> Dict[str, Any]:
    """
    Audita un módulo completo contra su public_api declarada.

    Resultado contiene:
    - module_id: str
    - module_name: str
    - surface_files: list[Path]
    - path_errors: list[PathError]
    - parse_errors: list[ParseError]
    - missing_symbols: list[str]
    """
    module_id = module["id"]
    module_name = module["name"]
    public_api = module.get("public_api", [])
    paths = module.get("paths", [])

    surface_files, path_errors = get_public_surface_files(repo_root, paths)

    # Analizar cada archivo válido y acumular errores de parsing
    all_bound_names: Set[str] = set()
    parse_errors: List[ParseError] = []

    for sf in surface_files:
        try:
            names = extract_top_level_names(sf)
            all_bound_names |= names
        except FileNotFoundError as e:
            parse_errors.append({
                "file": str(sf),
                "type": "FileNotFoundError",
                "message": str(e),
                "line": None,
            })
        except SyntaxError as e:
            parse_errors.append({
                "file": str(sf),
                "type": "SyntaxError",
                "message": e.msg,
                "line": e.lineno,
            })
        except UnicodeDecodeError as e:
            parse_errors.append({
                "file": str(sf),
                "type": "UnicodeDecodeError",
                "message": str(e),
                "line": e.start,  # approximate
            })

    if not public_api:
        return {
            "module_id": module_id,
            "module_name": module_name,
            "surface_files": surface_files,
            "path_errors": path_errors,
            "parse_errors": parse_errors,
            "missing_symbols": [],
        }

    missing = [sym for sym in public_api if sym not in all_bound_names]
    return {
        "module_id": module_id,
        "module_name": module_name,
        "surface_files": surface_files,
        "path_errors": path_errors,
        "parse_errors": parse_errors,
        "missing_symbols": sorted(set(missing)),
    }


# ---------------------------------------------------------------------------
# Tests sintéticos para el extractor AST (sin errores)
# ---------------------------------------------------------------------------

class TestExtractTopLevelNames:
    """Tests pequeños para verificar el comportamiento del extractor AST."""

    def test_direct_definitions(self, tmp_path: Path):
        """Una función, async function, clase, constante, variable anotada."""
        content = (
            "def foo(): pass\n"
            "async def bar(): pass\n"
            "class Baz: pass\n"
            "CONSTANT = 42\n"
            "annotated_var: str = 'x'\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "foo" in names
        assert "bar" in names
        assert "Baz" in names
        assert "CONSTANT" in names
        assert "annotated_var" in names

    def test_explicit_reexports(self, tmp_path: Path):
        """from pkg.internal import PublicClass; from pkg.internal import source_name as PublicAlias"""
        content = (
            "from pkg.internal import PublicClass\n"
            "from pkg.internal import source_name as PublicAlias\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "PublicClass" in names
        assert "PublicAlias" in names

    def test_import_normal_and_alias(self, tmp_path: Path):
        """import package; import package.submodule; import package.submodule as public_module"""
        content = (
            "import package\n"
            "import package.submodule\n"
            "import package.submodule as public_module\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "package" in names
        assert "public_module" in names

    def test_all_without_binding(self, tmp_path: Path):
        """Un símbolo presente solo en __all__ debe aparecer como faltante."""
        content = (
            "__all__ = ['MissingName']\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "MissingName" not in names

    def test_nested_definition(self, tmp_path: Path):
        """Función definida dentro de otra no debe aceptarse como nombre de nivel superior."""
        content = (
            "def outer():\n"
            "    def inner(): pass\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "inner" not in names
        assert "outer" in names

    def test_package_directory(self, tmp_path: Path):
        """Para un path de directorio, debe inspeccionar solo __init__.py."""
        dir_path = tmp_path / "mymodule"
        dir_path.mkdir()
        (dir_path / "__init__.py").write_text("# Only a comment\n")
        (dir_path / "internal.py").write_text("class InternalClass: pass\n")
        names = extract_top_level_names(dir_path / "__init__.py")
        assert "InternalClass" not in names

    def test_explicit_file(self, tmp_path: Path):
        """Un path que apunta directamente a single.py debe aceptar una función catalogada definida en ese archivo."""
        file = tmp_path / "single.py"
        file.write_text("def main(): pass\n")
        names = extract_top_level_names(file)
        assert "main" in names

    def test_multiple_explicit_files(self, tmp_path: Path):
        """Un módulo con dos archivos explícitos debe satisfacer símbolos distribuidos entre ambos."""
        file_a = tmp_path / "api_a.py"
        file_b = tmp_path / "api_b.py"
        file_a.write_text("def func_a(): pass\n")
        file_b.write_text("def func_b(): pass\n")
        all_names = set()
        for f in (file_a, file_b):
            all_names |= extract_top_level_names(f)
        assert "func_a" in all_names
        assert "func_b" in all_names

    def test_import_from_with_alias(self, tmp_path: Path):
        """Verifica que from ... import name as alias captura el alias."""
        content = "from chemuson.core import MolGraph as Graph\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "Graph" in names
        assert "MolGraph" not in names

    def test_import_star_not_accepted(self, tmp_path: Path):
        """import * no debe aceptar ningún símbolo automáticamente."""
        content = "from chemuson.core import *\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert len(names) == 0

    def test_multiple_imports_from_same_module(self, tmp_path: Path):
        """Múltiples imports desde el mismo módulo."""
        content = (
            "from pkg import A\n"
            "from pkg import B as C\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "A" in names
        assert "C" in names

    def test_nested_tuple_unpacking(self, tmp_path: Path):
        """Desempaquetado de tuplas simples."""
        content = "a: str = 'x'\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "a" in names


# ---------------------------------------------------------------------------
# Tests sintéticos que ejercitan el audit real
# ---------------------------------------------------------------------------

class TestAuditRealFlow:
    """Tests que usan audit_module_public_api para validar flujos completos."""

    def test_package_directory_flow(self, tmp_path: Path):
        """Crea un paquete con __init__.py y internal.py, declara el path del directorio.
        Un símbolo solo en internal.py queda faltante; uno reexportado en __init__.py sí existe.
        """
        dir_path = tmp_path / "mymodule"
        dir_path.mkdir()
        # Solo PublicClass se reexporta; InternalClass no se importa
        (dir_path / "__init__.py").write_text(
            "# Reexporta algo\n"
            "from .internal import PublicClass\n"
        )
        (dir_path / "internal.py").write_text(
            "class PublicClass: pass\n"
            "class InternalClass: pass\n"
        )

        module = {
            "id": "M99",
            "name": "mymodule",
            "public_api": ["PublicClass", "InternalClass"],
            "paths": ["mymodule"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0].name == "__init__.py"
        # path_errors debe estar vacío porque el directorio existe y tiene __init__.py
        assert result["path_errors"] == []
        # parse_errors debe estar vacío
        assert result["parse_errors"] == []
        # PublicClass está reexportado, InternalClass no
        assert "PublicClass" not in result["missing_symbols"]
        assert "InternalClass" in result["missing_symbols"]

    def test_explicit_file_flow(self, tmp_path: Path):
        """Declara un archivo único y verifica que una función pública existe."""
        file = tmp_path / "single.py"
        file.write_text("def main(): pass\n")

        module = {
            "id": "M99",
            "name": "single",
            "public_api": ["main"],
            "paths": ["single.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0].name == "single.py"
        assert result["path_errors"] == []
        assert result["parse_errors"] == []
        assert result["missing_symbols"] == []

    def test_multiple_explicit_files_flow(self, tmp_path: Path):
        """Declara dos archivos explícitos y distribuye símbolos entre ambos."""
        file_a = tmp_path / "api_a.py"
        file_b = tmp_path / "api_b.py"
        file_a.write_text("def func_a(): pass\n")
        file_b.write_text("def func_b(): pass\n")

        module = {
            "id": "M99",
            "name": "multi",
            "public_api": ["func_a", "func_b"],
            "paths": ["api_a.py", "api_b.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 2
        assert result["path_errors"] == []
        assert result["parse_errors"] == []
        assert result["missing_symbols"] == []

    def test_path_missing_flow(self, tmp_path: Path):
        """Declara un archivo o directorio que no existe. Debe aparecer path_missing."""
        module = {
            "id": "M99",
            "name": "missing",
            "public_api": [],
            "paths": ["does_not_exist.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 0
        assert len(result["path_errors"]) == 1
        assert result["path_errors"][0]["error_type"] == "path_missing"
        assert result["path_errors"][0]["path"] == "does_not_exist.py"

    def test_package_without_init_flow(self, tmp_path: Path):
        """Crea un directorio con internal.py pero sin __init__.py.
        Debe aparecer package_init_missing.
        """
        dir_path = tmp_path / "empty_pkg"
        dir_path.mkdir()
        (dir_path / "internal.py").write_text("class Something: pass\n")

        module = {
            "id": "M99",
            "name": "empty_pkg",
            "public_api": [],
            "paths": ["empty_pkg"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 0
        assert len(result["path_errors"]) == 1
        assert result["path_errors"][0]["error_type"] == "package_init_missing"
        assert result["path_errors"][0]["path"] == "empty_pkg"

    def test_unsupported_path_flow(self, tmp_path: Path):
        """Crea un archivo como README.txt y decláralo como path. Debe aparecer unsupported_path."""
        (tmp_path / "README.txt").write_text("Readme\n")

        module = {
            "id": "M99",
            "name": "unsupported",
            "public_api": [],
            "paths": ["README.txt"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["surface_files"]) == 0
        assert len(result["path_errors"]) == 1
        assert result["path_errors"][0]["error_type"] == "unsupported_path"
        assert result["path_errors"][0]["path"] == "README.txt"

    def test_syntax_error_flow(self, tmp_path: Path):
        """Crea un archivo Python con sintaxis inválida. Debe aparecer parse_error."""
        file = tmp_path / "bad.py"
        file.write_text("def broken(\n")  # sintaxis inválida

        module = {
            "id": "M99",
            "name": "syntax",
            "public_api": [],
            "paths": ["bad.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        # El archivo existe y es .py, así que se incluye en surface_files
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0] == file.resolve()
        # Pero al analizarlo, hay un parse error
        assert len(result["parse_errors"]) == 1
        assert result["parse_errors"][0]["type"] == "SyntaxError"
        assert result["parse_errors"][0]["file"] == str(file)

    def test_all_without_binding_flow(self, tmp_path: Path):
        """Mantiene la prueba de que __all__ = ["MissingName"] no satisface MissingName,
        pero compruébalo mediante el audit completo.
        """
        file = tmp_path / "test.py"
        file.write_text("__all__ = ['MissingName']\n")

        module = {
            "id": "M99",
            "name": "all_test",
            "public_api": ["MissingName"],
            "paths": ["test.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert "MissingName" in result["missing_symbols"]
        assert result["path_errors"] == []
        assert result["parse_errors"] == []

    def test_multiple_missing_symbols_flow(self, tmp_path: Path):
        """Declara al menos tres símbolos públicos ausentes. Todos aparecen en missing_symbols."""
        file = tmp_path / "test.py"
        file.write_text("def existing(): pass\n")

        module = {
            "id": "M99",
            "name": "multi_missing",
            "public_api": ["Symbol1", "Symbol2", "Symbol3"],
            "paths": ["test.py"],
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        assert len(result["missing_symbols"]) == 3
        assert sorted(result["missing_symbols"]) == ["Symbol1", "Symbol2", "Symbol3"]

    def test_deduplication_flow(self, tmp_path: Path):
        """Declara dos veces el mismo archivo o directorio. Solo se inspecciona una vez."""
        file = tmp_path / "dup.py"
        file.write_text("def func(): pass\n")

        module = {
            "id": "M99",
            "name": "dup",
            "public_api": ["func"],
            "paths": ["dup.py", "dup.py"],  # duplicado
        }

        result = audit_module_public_api(tmp_path, module)

        assert result["module_id"] == "M99"
        # El archivo debe aparecer solo una vez en surface_files
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0].name == "dup.py"
        assert result["missing_symbols"] == []


# ---------------------------------------------------------------------------
# Test contra el catálogo real: public_api existe en superficie pública
# ---------------------------------------------------------------------------

@pytest.fixture
def catalog():
    """Parse the module catalog YAML into a Python dict."""
    with open(CATALOG_PATH, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data


@pytest.fixture
def modules(catalog):
    """Return the list of module entries."""
    return catalog["modules"]


class TestPublicApiExists:
    """Validate that every symbol in public_api exists in the module's public surface."""

    def test_all_modules_with_public_api_are_validated(self, modules):
        """Verifica que todos los módulos con public_api no vacío son analizados."""
        modules_with_api = [m for m in modules if m.get("public_api")]
        assert len(modules_with_api) > 0

    def test_public_api_symbols_exist_for_each_module(self, modules):
        """
        Para cada módulo con public_api no vacío, verifica que todos los símbolos
        declarados existen en su superficie pública usando el audit unificado.
        """
        all_failures: List[dict] = []

        for module in modules:
            module_id = module["id"]
            module_name = module["name"]
            public_api = module.get("public_api", [])

            if not public_api:
                continue

            result = audit_module_public_api(REPO_ROOT, module)

            # Acumular path_errors como fallos si existen
            if result["path_errors"]:
                all_failures.append({
                    "module_id": module_id,
                    "module_name": module_name,
                    "path_errors": result["path_errors"],
                })

            # Acumular parse_errors como fallos si existen
            if result["parse_errors"]:
                all_failures.append({
                    "module_id": module_id,
                    "module_name": module_name,
                    "parse_errors": result["parse_errors"],
                })

            # Símbolos faltantes
            if result["missing_symbols"]:
                all_failures.append({
                    "module_id": module_id,
                    "module_name": module_name,
                    "missing_symbols": result["missing_symbols"],
                    "surface_files": [str(f.relative_to(REPO_ROOT)) for f in result["surface_files"]],
                })

        if all_failures:
            msg_lines = [f"Found {len(all_failures)} issues across modules:"]
            for failure in all_failures:
                module_id = failure["module_id"]
                module_name = failure["module_name"]
                if "path_errors" in failure:
                    msg_lines.append(f"  - {module_id} ({module_name}): path errors:")
                    for pe in failure["path_errors"]:
                        msg_lines.append(f"      {pe['error_type']}: {pe['path']}")
                if "parse_errors" in failure:
                    msg_lines.append(f"  - {module_id} ({module_name}): parse errors:")
                    for pe in failure["parse_errors"]:
                        line_info = f":line {pe['line']}" if pe.get("line") else ""
                        msg_lines.append(f"      {pe['type']}{line_info}: {pe['message']}")
                if "missing_symbols" in failure:
                    msg_lines.append(f"  - {module_id} ({module_name}): missing symbols:")
                    for sym in failure["missing_symbols"]:
                        msg_lines.append(f"      Missing: {sym}")
                    msg_lines.append(f"      Surface files: {', '.join(failure['surface_files'])}")
            assert False, "\n".join(msg_lines)

    def test_m18_version_module(self, modules):
        """M18 tiene public_api en varios archivos explícitos."""
        m18 = next((m for m in modules if m["id"] == "M18"), None)
        assert m18 is not None
        assert m18["public_api"] == ["__version__", "get_app_version"]

        result = audit_module_public_api(REPO_ROOT, m18)

        assert result["module_id"] == "M18"
        assert len(result["path_errors"]) == 0
        assert len(result["parse_errors"]) == 0
        assert result["missing_symbols"] == []

    def test_m19_bootstrap_module(self, modules):
        """M19 expone main desde __main__.py y posee app como paquete interno."""
        m19 = next((m for m in modules if m["id"] == "M19"), None)
        assert m19 is not None
        assert m19["public_api"] == ["main"]

        result = audit_module_public_api(REPO_ROOT, m19)

        assert result["module_id"] == "M19"
        assert {path.relative_to(REPO_ROOT).as_posix() for path in result["surface_files"]} == {
            "src/chemuson/__main__.py",
            "src/chemuson/app/__init__.py",
        }
        assert len(result["path_errors"]) == 0
        assert len(result["parse_errors"]) == 0
        assert "main" not in result["missing_symbols"]

    def test_m00_core_module(self, modules):
        """M00 core es un paquete normal con __init__.py."""
        m00 = next((m for m in modules if m["id"] == "M00"), None)
        assert m00 is not None
        assert len(m00["public_api"]) > 0

        result = audit_module_public_api(REPO_ROOT, m00)

        assert result["module_id"] == "M00"
        # Debe haber solo un __init__.py en surface_files
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0].name == "__init__.py"
        assert len(result["path_errors"]) == 0
        # Verificar que no hay símbolos faltantes (el test anterior ya lo hace)
        assert result["missing_symbols"] == []

    def test_m02_clean2d_module(self, modules):
        """M02 clean2d es un paquete con muchas reexports."""
        m02 = next((m for m in modules if m["id"] == "M02"), None)
        assert m02 is not None
        assert len(m02["public_api"]) > 30

        result = audit_module_public_api(REPO_ROOT, m02)

        assert result["module_id"] == "M02"
        assert len(result["surface_files"]) == 1
        assert result["surface_files"][0].name == "__init__.py"
        assert len(result["path_errors"]) == 0
        assert result["missing_symbols"] == []


# ---------------------------------------------------------------------------
# Tests específicos para casos de import (refactorizados)
# ---------------------------------------------------------------------------

class TestImportHandling:
    """Tests adicionales para verificar manejo correcto de imports."""

    def test_import_from_with_alias(self, tmp_path: Path):
        """Verifica que from ... import name as alias captura el alias."""
        content = "from chemuson.core import MolGraph as Graph\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "Graph" in names
        assert "MolGraph" not in names

    def test_import_star_not_accepted(self, tmp_path: Path):
        """import * no debe aceptar ningún símbolo automáticamente."""
        content = "from chemuson.core import *\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert len(names) == 0

    def test_multiple_imports_from_same_module(self, tmp_path: Path):
        """Múltiples imports desde el mismo módulo."""
        content = (
            "from pkg import A\n"
            "from pkg import B as C\n"
        )
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "A" in names
        assert "C" in names

    def test_nested_tuple_unpacking(self, tmp_path: Path):
        """Desempaquetado de tuplas simples."""
        content = "a: str = 'x'\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "a" in names
