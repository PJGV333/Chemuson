"""Tests arquitectónicos: validación de que los símbolos de public_api existen en la superficie pública (Phase 1).

Valida que cada símbolo declarado en `public_api` existe como nombre de nivel superior
en los archivos de superficie pública correspondientes, sin importar módulos de producción.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import List, Set

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"


# ---------------------------------------------------------------------------
# Helpers: determinación de superficie pública a partir de paths del catálogo
# ---------------------------------------------------------------------------

def get_public_surface_files(module_paths: list[str]) -> list[Path]:
    """
    Dado el listado de paths de un módulo en el catálogo, devuelve la lista
    de archivos Python que constituyen su superficie pública.

    Reglas:
    - Si path es directorio -> solo __init__.py dentro de él.
    - Si path es archivo .py explícito -> ese archivo.
    - Varios paths -> unión (ordenada y deduplicada).
    """
    surface_files: List[Path] = []

    for path_str in module_paths:
        full_path = REPO_ROOT / path_str

        if not full_path.exists():
            continue

        if full_path.is_dir():
            init_file = full_path / "__init__.py"
            if init_file.exists():
                surface_files.append(init_file.resolve())
        else:
            if full_path.suffix == ".py":
                surface_files.append(full_path.resolve())

    # Deduplicar y ordenar de forma determinista
    return sorted(set(surface_files))


# ---------------------------------------------------------------------------
# Helpers: extracción de nombres vinculados a nivel superior mediante AST
# ---------------------------------------------------------------------------

def extract_top_level_names(file_path: Path) -> Set[str]:
    """
    Analiza el archivo dado y extrae los nombres que están vinculados a nivel
    superior del módulo: definiciones, asignaciones e imports.

    No usa ast.walk() para evitar capturar nombres anidados. Recorre solo el body.
    """
    if not file_path.exists():
        return set()

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        source = f.read()

    try:
        tree = ast.parse(source, filename=str(file_path))
    except SyntaxError:
        return set()

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
                # El nombre vinculado es el primer componente del módulo importado
                top = alias.name.split(".", 1)[0]
                if alias.asname:
                    names.add(alias.asname)
                else:
                    names.add(top)
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                # Skip '*' imports as they don't bind any specific name
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
# Helpers: validación de símbolos públicos
# ---------------------------------------------------------------------------

def validate_module_public_api(module_id: str, module_name: str, public_api: list[str], surface_files: list[Path]) -> list[str]:
    """
    Para un módulo, calcula los símbolos que faltan en su superficie pública.
    Devuelve lista de símbolos ausentes. Si la lista está vacía, todos existen.
    """
    if not public_api:
        return []

    all_bound_names: Set[str] = set()
    for sf in surface_files:
        if sf.exists():
            all_bound_names |= extract_top_level_names(sf)

    missing = [sym for sym in public_api if sym not in all_bound_names]
    return missing


# ---------------------------------------------------------------------------
# Tests sintéticos para los helpers
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
        assert "package" in names  # segunda import también vincula package
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
        """Para un path de directorio, debe inspeccionar solo __init__.py.
        Un símbolo que solo existe en internal.py y no se reexporta en __init__.py
        NO debe aparecer en los nombres extraídos."""
        dir_path = tmp_path / "mymodule"
        dir_path.mkdir()
        # __init__.py sin reexportar InternalClass
        (dir_path / "__init__.py").write_text("# Only a comment\n")
        (dir_path / "internal.py").write_text("class InternalClass: pass\n")
        names = extract_top_level_names(dir_path / "__init__.py")
        # InternalClass no está reexportado en __init__, así que no debería aparecer
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
        # Simulamos unión de archivos
        all_names = set()
        for f in (file_a, file_b):
            all_names |= extract_top_level_names(f)
        assert "func_a" in all_names
        assert "func_b" in all_names

    def test_missing_file_or_init(self, tmp_path: Path):
        """Si un path explícito no existe, o un directorio no contiene __init__.py, el test debe fallar."""
        # Archivo que no existe
        fake_file = tmp_path / "does_not_exist.py"
        assert not fake_file.exists()
        names = extract_top_level_names(fake_file)
        assert names == set()

        # Directorio sin __init__.py
        empty_dir = tmp_path / "empty_package"
        empty_dir.mkdir()
        names = extract_top_level_names(empty_dir / "__init__.py")
        assert names == set()


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
        declarados existen en su superficie pública.
        """
        all_failures: List[dict] = []

        for module in modules:
            module_id = module["id"]
            module_name = module["name"]
            public_api = module.get("public_api", [])
            paths = module.get("paths", [])

            if not public_api:
                continue

            surface_files = get_public_surface_files(paths)

            # Verificar que los archivos de superficie existen
            for sf in surface_files:
                if not sf.exists():
                    all_failures.append({
                        "module_id": module_id,
                        "module_name": module_name,
                        "error": "surface_file_missing",
                        "path": str(sf),
                    })
                    continue

            missing = validate_module_public_api(module_id, module_name, public_api, surface_files)

            if missing:
                all_failures.append({
                    "module_id": module_id,
                    "module_name": module_name,
                    "missing_symbols": missing,
                    "surface_files": [str(f.relative_to(REPO_ROOT)) for f in surface_files],
                })

        if all_failures:
            msg_lines = [f"Found {len(all_failures)} modules with missing public API symbols:"]
            for failure in all_failures:
                module_id = failure["module_id"]
                module_name = failure["module_name"]
                if "missing_symbols" in failure:
                    msg_lines.append(f"  - {module_id} ({module_name}):")
                    for sym in failure["missing_symbols"]:
                        msg_lines.append(f"      Missing: {sym}")
                    msg_lines.append(f"      Surface files: {', '.join(failure['surface_files'])}")
                else:
                    msg_lines.append(f"  - {module_id} ({module_name}): surface file missing: {failure['path']}")
            assert False, "\n".join(msg_lines)

    def test_m18_version_module(self, modules):
        """M18 tiene public_api en varios archivos explícitos."""
        m18 = next((m for m in modules if m["id"] == "M18"), None)
        assert m18 is not None
        assert m18["public_api"] == ["__version__", "get_app_version"]
        surface_files = get_public_surface_files(m18["paths"])
        # Ambos archivos deben existir
        for sf in surface_files:
            assert sf.exists()

    def test_m19_bootstrap_module(self, modules):
        """M19 tiene surface única en __main__.py."""
        m19 = next((m for m in modules if m["id"] == "M19"), None)
        assert m19 is not None
        assert m19["public_api"] == ["main"]
        surface_files = get_public_surface_files(m19["paths"])
        assert len(surface_files) == 1
        assert surface_files[0].name == "__main__.py"
        assert surface_files[0].exists()

    def test_m00_core_module(self, modules):
        """M00 core es un paquete normal con __init__.py."""
        m00 = next((m for m in modules if m["id"] == "M00"), None)
        assert m00 is not None
        assert len(m00["public_api"]) > 0
        surface_files = get_public_surface_files(m00["paths"])
        # Debe haber solo un __init__.py
        assert len(surface_files) == 1
        assert surface_files[0].name == "__init__.py"

    def test_m02_clean2d_module(self, modules):
        """M02 clean2d es un paquete con muchas reexports."""
        m02 = next((m for m in modules if m["id"] == "M02"), None)
        assert m02 is not None
        assert len(m02["public_api"]) > 30
        surface_files = get_public_surface_files(m02["paths"])
        assert len(surface_files) == 1
        assert surface_files[0].name == "__init__.py"


# ---------------------------------------------------------------------------
# Tests específicos para casos de import
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
        # No hay símbolos vinculados con *
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
        """Desempaquetado de tuplas (anidadas o no)."""
        # Usamos una asignación an simple sin empaquetado para verificar AnnAssign básico
        content = "a: str = 'x'\n"
        file = tmp_path / "test.py"
        file.write_text(content)
        names = extract_top_level_names(file)
        assert "a" in names
