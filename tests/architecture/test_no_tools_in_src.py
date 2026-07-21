"""Tests arquitectónicos: verificación de que src/chemuson/ no importa el paquete tools.

Verifica que ningún archivo Python bajo `src/chemuson/` importe desde `tools` o `chemuson.tools`.
Esta regla no tiene excepciones, ni siquiera para TYPE_CHECKING.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent.parent.parent / "src" / "chemuson"


@dataclass(frozen=True)
class ToolsImportViolation:
    """Representa una violación detectada de import prohibido a tools."""
    relative_file: str
    line: int
    col: int
    statement: str
    forbidden_module: str
    node_type: str  # "import" o "from"
    rule: str = "production_must_not_import_tools"

    def __lt__(self, other: "ToolsImportViolation") -> bool:
        return (
            self.relative_file,
            self.line,
            self.col,
            self.forbidden_module,
        ) < (
            other.relative_file,
            other.line,
            other.col,
            other.forbidden_module,
        )


def is_forbidden_tools_module(module_name: str) -> bool:
    """
    Determina si un nombre de módulo importado corresponde a `tools` o `chemuson.tools`.

    Reglas:
    - module_name == "tools" es prohibido.
    - module_name empieza con "tools." es prohibido.
    - module_name == "chemuson.tools" es prohibido.
    - module_name empieza con "chemuson.tools." es prohibido.

    No se usan coincidencias por subcadena (ej: "tools" en nombre) para evitar falsos positivos.
    """
    if module_name == "tools":
        return True
    if module_name.startswith("tools."):
        return True
    if module_name == "chemuson.tools":
        return True
    if module_name.startswith("chemuson.tools."):
        return True
    return False


def get_module_name_for_file(src_root: Path, file_path: Path) -> str:
    """
    Calcula el nombre de módulo completo para un archivo dentro del árbol de src/chemuson/.

    Ejemplos:
    - src/chemuson/core/model.py -> chemuson.core.model
    - src/chemuson/gui/__init__.py -> chemuson.gui
    - src/chemuson/__init__.py -> chemuson
    """
    if not file_path.is_relative_to(src_root):
        raise ValueError(f"File not under src root: {file_path}")

    rel_path = file_path.relative_to(src_root)
    parts = list(rel_path.parts)

    # Eliminar __init__.py del nombre final
    if parts[-1] == "__init__.py":
        parts = parts[:-1]
    elif parts[-1].endswith(".py"):
        parts[-1] = parts[-1][: -len(".py")]

    if not parts:
        return "chemuson"

    return "chemuson." + ".".join(parts)


def resolve_import_from_base(
    base_module: str,
    level: int,
    module: str | None,
) -> str | None:
    """
    Resuelve el módulo completo a partir de un import relativo.

    level: ImportFrom.level (1 = from . import, 2 = from .. import, etc.)
    base_module: módulo del archivo que realiza el import
    module: módulo especificado en el import (puede ser None)

    Devuelve el módulo completo resuelto o None si no se puede resolver.
    """
    if level == 0:
        # Import absoluto, el módulo ya es completo
        return module

    parts = base_module.split(".")

    if level > len(parts):
        # Nivel de import relativo supera la profundidad del paquete
        return None

    # Calcular la parte base después de aplicar el nivel relativo
    # Para un archivo en chemuson.gui.controllers, level=1 -> base es chemuson.gui.controllers
    # Para level=2 -> base es chemuson.gui
    if parts[-1] == "__init__":
        # Si el archivo es __init__.py, el paquete es la parte padre
        base_parts = parts[: -level + 1] if level > 1 else parts
    else:
        # El módulo base es el contenedor del archivo
        base_parts = parts[: -level]

    if not base_parts:
        base_parts = ["chemuson"]

    if module:
        return ".".join(base_parts + [module])
    else:
        return ".".join(base_parts)


def scan_file_for_tools_imports(
    file_path: Path,
    src_root: Path,
) -> tuple[list[ToolsImportViolation], list[dict], int]:
    """
    Escanea un archivo para detectar imports prohibidos de `tools` o `chemuson.tools`.

    Devuelve:
    - lista de violaciones encontradas
    - lista de errores de parsing (si hay)
    """
    violations: list[ToolsImportViolation] = []
    parse_errors: list[dict] = []

    try:
        with open(file_path, "r", encoding="utf-8", errors="replace") as f:
            source = f.read()
    except Exception as e:
        parse_errors.append({
            "file": str(file_path),
            "type": "FileReadError",
            "message": str(e),
            "line": None,
            "col": None,
        })
        return violations, parse_errors, 0

    try:
        tree = ast.parse(source, filename=str(file_path))
    except SyntaxError as e:
        parse_errors.append({
            "file": str(file_path),
            "type": "SyntaxError",
            "message": e.msg,
            "line": e.lineno,
            "col": e.offset,
        })
        return violations, parse_errors, 0

    # Calcular el módulo base del archivo
    try:
        base_module = get_module_name_for_file(src_root, file_path)
    except ValueError as e:
        parse_errors.append({
            "file": str(file_path),
            "type": "ValueError",
            "message": str(e),
            "line": None,
            "col": None,
        })
        return violations, parse_errors, 0

    # Recorrer todos los nodos del AST para encontrar imports en cualquier profundidad
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                # Para `import tools` o `import tools.release`
                if is_forbidden_tools_module(alias.name):
                    statement = ast.get_source_segment(source, node)
                    if statement is None:
                        statement = ast.unparse(node)
                    statement = " ".join(statement.split())
                    violations.append(ToolsImportViolation(
                        relative_file=str(file_path.relative_to(src_root.parent)),
                        line=node.lineno,
                        col=node.col_offset,
                        statement=statement,
                        forbidden_module=alias.name,
                        node_type="import",
                    ))

        elif isinstance(node, ast.ImportFrom):
            # Resolver el módulo completo
            resolved = resolve_import_from_base(base_module, node.level, node.module)

            # Primero, verificar si el módulo resuelto es prohibido
            if resolved and is_forbidden_tools_module(resolved):
                for alias in node.names:
                    statement = ast.get_source_segment(source, node)
                    if statement is None:
                        statement = ast.unparse(node)
                    statement = " ".join(statement.split())
                    violations.append(ToolsImportViolation(
                        relative_file=str(file_path.relative_to(src_root.parent)),
                        line=node.lineno,
                        col=node.col_offset,
                        statement=statement,
                        forbidden_module=resolved,
                        node_type="from",
                    ))

            # Segundo, verificar patrones especiales donde el nombre importado es "tools"
            # y podría ser un submodule bajo chemuson.tools

            # Caso: from <pkg> import tools (import de nivel 0 sin alias del módulo)
            # Ejemplo: from chemuson import tools -> candidato: chemuson.tools
            if node.module and node.level == 0:
                for alias in node.names:
                    if alias.name == "tools":
                        candidate = f"{node.module}.tools"
                        if is_forbidden_tools_module(candidate):
                            statement = ast.get_source_segment(source, node)
                            if statement is None:
                                statement = ast.unparse(node)
                            statement = " ".join(statement.split())
                            violations.append(ToolsImportViolation(
                                relative_file=str(file_path.relative_to(src_root.parent)),
                                line=node.lineno,
                                col=node.col_offset,
                                statement=statement,
                                forbidden_module=candidate,
                                node_type="from",
                            ))

            # Caso: from . import tools (relative sin módulo explícito)
            # Ejemplo: from . import tools -> candidato: <base_package>.tools
            if node.module is None and node.level > 0:
                for alias in node.names:
                    if alias.name == "tools":
                        parts = base_module.split(".")
                        if node.level > len(parts):
                            continue
                        if parts[-1] == "__init__":
                            base_parts = parts[: -node.level + 1] if node.level > 1 else parts
                        else:
                            base_parts = parts[: -node.level]
                        if not base_parts:
                            base_parts = ["chemuson"]
                        candidate = ".".join(base_parts + ["tools"])
                        if is_forbidden_tools_module(candidate):
                            statement = ast.get_source_segment(source, node)
                            if statement is None:
                                statement = ast.unparse(node)
                            statement = " ".join(statement.split())
                            violations.append(ToolsImportViolation(
                                relative_file=str(file_path.relative_to(src_root.parent)),
                                line=node.lineno,
                                col=node.col_offset,
                                statement=statement,
                                forbidden_module=candidate,
                                node_type="from",
                            ))

            # Tercero, verificar if .tools import (relative con módulo)
            # Ejemplo: from .tools import helper -> module="tools", level=1
            # Esto ya es manejado por resolve_import_from_base?
            # Para from .tools import helper, node.module="tools", level=1
            # La resolución debería dar "chemuson.tools" si base_module="chemuson"
            # Pero si base_module="chemuson.gui", nivel=1, módulo="tools" -> chemuson.tools
            # Eso ya está cubierto por el primer bloque (resolved).


    return violations, parse_errors, 1


def scan_src_for_tools_imports(src_root: Path) -> tuple[list[ToolsImportViolation], list[dict], int]:
    """
    Escanea recursivamente todos los archivos .py bajo src_root.

    Devuelve:
    - lista de todas las violaciones encontradas
    - lista de todos los errores de parsing
    - número de archivos analizados
    """
    all_violations: list[ToolsImportViolation] = []
    all_parse_errors: list[dict] = []
    files_analyzed = 0

    for file_path in src_root.rglob("*.py"):
        if file_path.name == "__pycache__":
            continue
        files_analyzed += 1
        violations, parse_errors, _ = scan_file_for_tools_imports(file_path, src_root)
        all_violations.extend(violations)
        all_parse_errors.extend(parse_errors)

    # Ordenar resultados de forma determinista
    all_violations.sort()
    all_parse_errors.sort(key=lambda e: (e["file"], e.get("line", 0), e.get("col", 0)))

    return all_violations, all_parse_errors, files_analyzed


class TestNoToolsImports:
    """Tests arquitectónicos para verificar que src/chemuson/ no importa herramientas."""

    def test_src_directory_exists(self):
        """Debe existir el directorio src/chemuson/."""
        assert SRC_ROOT.exists(), f"Directorio no encontrado: {SRC_ROOT}"
        assert SRC_ROOT.is_dir(), f"La ruta no es un directorio: {SRC_ROOT}"

    def test_found_files_to_analyze(self):
        """Debe encontrar al menos un archivo .py para analizar."""
        file_count = sum(1 for _ in SRC_ROOT.rglob("*.py") if _.name != "__pycache__")
        assert file_count > 0, f"No se encontraron archivos .py bajo {SRC_ROOT}"

    def test_no_tools_imports_in_production(self):
        """
        El test principal: ningún archivo en src/chemuson/ debe importar tools o chemuson.tools.

        Se acumulan todos los errores en una sola ejecución.
        """
        violations, parse_errors, files_count = scan_src_for_tools_imports(SRC_ROOT)

        # Primero reportar cualquier error de parsing (ocultaría violaciones)
        if parse_errors:
            msg_lines = [f"Found {len(parse_errors)} parsing errors that prevent analysis:"]
            for err in parse_errors:
                line_info = f":line {err['line']}" if err.get("line") else ""
                col_info = f":col {err['col']}" if err.get("col") else ""
                msg_lines.append(f"  - {err['file']}{line_info}{col_info}: {err['type']} - {err['message']}")
            assert False, "\n".join(msg_lines)

        # Luego reportar violaciones
        if violations:
            msg_lines = [
                f"Found {len(violations)} forbidden tools imports across {files_count} analyzed files:",
            ]
            for v in violations:
                msg_lines.append(
                    f"  - {v.relative_file}:{v.line}:{v.col} | "
                    f"{v.node_type} | {v.forbidden_module} | {v.statement}"
                )
            assert False, "\n".join(msg_lines)


# ---------------------------------------------------------------------------
# Tests sintéticos para los helpers
# ---------------------------------------------------------------------------

class TestIsForbiddenToolsModule:
    """Tests unitarios para la función de detección de módulos prohibidos."""

    def test_tools_exact(self):
        assert is_forbidden_tools_module("tools") is True

    def test_tools_with_submodule(self):
        assert is_forbidden_tools_module("tools.release") is True
        assert is_forbidden_tools_module("tools.chemname_acceptance") is True

    def test_chemuson_tools_exact(self):
        assert is_forbidden_tools_module("chemuson.tools") is True

    def test_chemuson_tools_with_submodule(self):
        assert is_forbidden_tools_module("chemuson.tools.release") is True
        assert is_forbidden_tools_module("chemuson.tools.visual") is True

    def test_similar_names_not_forbidden(self):
        # Estos NO deben ser prohibidos (falsos positivos)
        assert is_forbidden_tools_module("tooling") is False
        assert is_forbidden_tools_module("toolshed") is False
        assert is_forbidden_tools_module("chemuson_tooling") is False
        assert is_forbidden_tools_module("chemuson.tools_helper") is False
        assert is_forbidden_tools_module("toolsmith") is False


class TestGetModuleNameForFile:
    """Tests unitarios para el cálculo del nombre de módulo."""

    def test_file_in_package(self, tmp_path: Path):
        src_root = tmp_path / "src" / "chemuson"
        file = src_root / "core" / "model.py"
        file.parent.mkdir(parents=True)
        file.write_text("# test")
        result = get_module_name_for_file(src_root, file)
        assert result == "chemuson.core.model"

    def test_init_file_in_package(self, tmp_path: Path):
        src_root = tmp_path / "src" / "chemuson"
        file = src_root / "gui" / "__init__.py"
        file.parent.mkdir(parents=True)
        file.write_text("# test")
        result = get_module_name_for_file(src_root, file)
        assert result == "chemuson.gui"

    def test_init_file_at_root(self, tmp_path: Path):
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "__init__.py"
        file.write_text("# test")
        result = get_module_name_for_file(src_root, file)
        assert result == "chemuson"


class TestResolveImportFromBase:
    """Tests unitarios para la resolución de imports relativos."""

    def test_absolute_import(self):
        assert resolve_import_from_base("chemuson.core", 0, "tools") == "tools"
        assert resolve_import_from_base("chemuson.core", 0, "chemuson.tools") == "chemuson.tools"

    def test_relative_level_1_from_package_file(self):
        # Archivo: chemuson/gui/controllers/file.py
        # from . import tools -> level=1, base should be chemuson.gui.controllers
        assert resolve_import_from_base("chemuson.gui.controllers", 1, "tools") == "chemuson.gui.tools"

    def test_relative_level_2_from_package_file(self):
        # Archivo: chemuson/gui/controllers/file.py
        # from .. import tools -> level=2, base should be chemuson.gui
        assert resolve_import_from_base("chemuson.gui.controllers", 2, "tools") == "chemuson.tools"

    def test_relative_level_1_from_init_file(self):
        # Archivo: chemuson/gui/__init__.py
        # from . import tools -> level=1, base should be chemuson.gui (after removing __init__)
        # But in get_module_name_for_file, __init__.py is removed, so base is chemuson.gui
        assert resolve_import_from_base("chemuson.gui", 1, "tools") == "chemuson.tools"

    def test_relative_level_2_from_subpackage_init(self):
        # Archivo: chemuson/gui/canvas/__init__.py
        # from .. import tools -> level=2, base should be chemuson
        assert resolve_import_from_base("chemuson.gui.canvas", 2, "tools") == "chemuson.tools"


# ---------------------------------------------------------------------------
# Tests sintéticos completos que usan el flujo real (scan_file_for_tools_imports)
# ---------------------------------------------------------------------------

class TestScanFileRealFlow:
    """Tests que usan scan_file_for_tools_imports para verificar flujos completos."""

    def test_import_tools(self, tmp_path: Path):
        """Debe producir una violación para `import tools`."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("import tools\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 1
        assert violations[0].forbidden_module == "tools"
        assert violations[0].node_type == "import"

    def test_submodule_and_alias(self, tmp_path: Path):
        """Debe detectar `import tools.release as release_tools`."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("import tools.release as release_tools\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 1
        assert violations[0].forbidden_module == "tools.release"

    def test_from_tools(self, tmp_path: Path):
        """Debe detectar `from tools import helper` y `from tools.release import build`."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("from tools import helper\nfrom tools.release import build\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 2
        forbidden_modules = {v.forbidden_module for v in violations}
        assert forbidden_modules == {"tools", "tools.release"}

    def test_chemuson_tools(self, tmp_path: Path):
        """Debe detectar `import chemuson.tools` y `from chemuson.tools.release import build`."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("import chemuson.tools\nfrom chemuson.tools.release import build\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 2
        forbidden_modules = {v.forbidden_module for v in violations}
        assert forbidden_modules == {"chemuson.tools", "chemuson.tools.release"}

    def test_from_chemuson_import_tools(self, tmp_path: Path):
        """Debe detectar `from chemuson import tools` y con alias."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("from chemuson import tools\nfrom chemuson import tools as repo_tools\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 2
        assert all(v.forbidden_module == "chemuson.tools" for v in violations)

    def test_relative_from_package_root(self, tmp_path: Path):
        """Debe detectar imports relativos desde la raíz del paquete."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text("from .tools import helper\nfrom . import tools\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        # El archivo test.py está directamente bajo src/chemuson, base_module = chemuson
        # from .tools import helper -> level=1, module=tools -> chemuson.tools
        assert len(violations) == 2
        assert all(v.forbidden_module == "chemuson.tools" for v in violations)

    def test_relative_from_subpackage(self, tmp_path: Path):
        """Debe detectar imports relativos desde subpaquetes."""
        src_root = tmp_path / "src" / "chemuson"
        subdir = src_root / "gui"
        subdir.mkdir(parents=True)
        file = subdir / "test.py"
        file.write_text("from ..tools import helper\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        # Base module: chemuson.gui (file is under gui/)
        # from ..tools import helper -> level=2, module=tools -> chemuson.tools
        assert len(violations) == 1
        assert violations[0].forbidden_module == "chemuson.tools"

    def test_import_in_function_and_type_checking(self, tmp_path: Path):
        """Confirma que ambos se detectan: imports dentro de funciones y TYPE_CHECKING."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text(
            "def function():\n"
            "    import tools\n"
            "if TYPE_CHECKING:\n"
            "    from chemuson.tools import ToolType\n"
        )

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 2
        forbidden_modules = {v.forbidden_module for v in violations}
        assert forbidden_modules == {"tools", "chemuson.tools"}

    def test_similar_names_allowed(self, tmp_path: Path):
        """Confirma cero violaciones para nombres parecidos permitidos."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text(
            "import tooling\n"
            "import toolshed\n"
            "from chemuson import tools_helper\n"
        )

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 0

    def test_canvas_tools_allowed(self, tmp_path: Path):
        """Confirma cero violaciones para imports como canvas_tools_bonding."""
        src_root = tmp_path / "src" / "chemuson"
        subdir = src_root / "gui" / "canvas"
        subdir.mkdir(parents=True)
        file = subdir / "test.py"
        file.write_text(
            "from .canvas_tools_bonding import CanvasToolsBondingMixin\n"
            "from chemuson.gui.canvas.canvas_tools_annotations import CanvasToolsAnnotationsMixin\n"
        )

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 0

    def test_strings_and_comments(self, tmp_path: Path):
        """Este contenido no debe producir violaciones."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text(
            "# import tools\n"
            "TEXT = 'from chemuson.tools import helper'\n"
        )

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 0

    def test_multiple_violations(self, tmp_path: Path):
        """Un archivo con varias violaciones debe devolverlas todas en orden determinista."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "test.py"
        file.write_text(
            "import tools\n"
            "from tools import helper\n"
            "import chemuson.tools\n"
            "from chemuson.tools import build\n"
        )

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 4
        # Verificar orden
        assert violations[0].line == 1
        assert violations[3].line == 4

    def test_syntax_error(self, tmp_path: Path):
        """Un archivo con sintaxis inválida debe producir un error de parsing explícito."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        file = src_root / "bad.py"
        file.write_text("def broken(\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        assert len(violations) == 0
        assert len(parse_errors) == 1
        assert parse_errors[0]["type"] == "SyntaxError"

    def test_clean_files(self, tmp_path: Path):
        """Un árbol sintético con varios .py sin imports prohibidos debe producir cero violaciones."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        (src_root / "a.py").write_text("def f(): pass\n")
        (src_root / "b.py").write_text("from .a import f\n")
        subdir = src_root / "sub"
        subdir.mkdir()
        (subdir / "c.py").write_text("from ..a import f\n")

        violations, parse_errors, files_count = scan_src_for_tools_imports(src_root)

        assert len(violations) == 0
        assert len(parse_errors) == 0
        assert files_count == 3

    def test_init_file_resolution(self, tmp_path: Path):
        """Comprueba que el paquete actual de src/chemuson/gui/__init__.py se calcula como chemuson.gui.
        from ..tools import something debe importar chemuson.tools."""
        src_root = tmp_path / "src" / "chemuson"
        src_root.mkdir(parents=True)
        gui_dir = src_root / "gui"
        gui_dir.mkdir(parents=True)
        file = gui_dir / "__init__.py"
        file.write_text("from ..tools import something\n")

        violations, parse_errors, _ = scan_file_for_tools_imports(file, src_root)

        # Base module: chemuson.gui (not chemuson.gui.__init__)
        # from ..tools import something -> level=2, module=tools -> chemuson.tools
        assert len(violations) == 1
        assert violations[0].forbidden_module == "chemuson.tools"
