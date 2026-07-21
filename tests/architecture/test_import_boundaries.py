import ast
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

import pytest
import yaml


@dataclass(frozen=True)
class ObservedImport:
    source_id: str
    target_id: str
    file: Path
    line: int
    statement: str
    type_checking_only: bool


class ImportVisitor(ast.NodeVisitor):
    def __init__(self, source_text: str, file_path: Path, current_module_name: str):
        self.source_text = source_text
        self.file_path = file_path
        self.current_module_name = current_module_name
        self.imports: List[Dict[str, Any]] = []
        self.in_type_checking = False
        self.context_stack: List[bool] = []

    def _is_type_checking_test(self, test_node: ast.expr) -> bool:
        if isinstance(test_node, ast.Name) and test_node.id == "TYPE_CHECKING":
            return True
        if (isinstance(test_node, ast.Attribute) and
                isinstance(test_node.value, ast.Name) and
                test_node.value.id == "typing" and
                test_node.attr == "TYPE_CHECKING"):
            return True
        return False

    def _visit_body(self, body: list[ast.stmt]) -> None:
        for stmt in body:
            self.visit(stmt)

    def visit_If(self, node: ast.If):
        is_tc = self._is_type_checking_test(node.test)
        prev_tc = self.in_type_checking
        self.context_stack.append(prev_tc)

        if is_tc:
            self.in_type_checking = True
            self._visit_body(node.body)
            self.in_type_checking = prev_tc
            self._visit_body(node.orelse)
        else:
            self._visit_body(node.body)
            self._visit_body(node.orelse)

        self.in_type_checking = self.context_stack.pop()

    def visit_Import(self, node: ast.Import):
        for alias in node.names:
            self.imports.append({
                "type": "absolute",
                "name": alias.name,
                "lineno": node.lineno,
                "node": node,
                "type_checking_only": self.in_type_checking
            })

    def visit_ImportFrom(self, node: ast.ImportFrom):
        self.imports.append({
            "type": "from",
            "module": node.module if node.module else "",
            "level": node.level,
            "names": [alias.name for alias in node.names],
            "lineno": node.lineno,
            "node": node,
            "type_checking_only": self.in_type_checking
        })


class ImportBoundaryAnalyzer:
    def __init__(self, repo_root: Path, catalog_path: Path):
        self.repo_root = repo_root
        self.catalog_path = catalog_path
        self.modules_cfg = self._load_catalog()
        self.file_to_module: Dict[Path, str] = {}
        self.importable_to_module_id: Dict[str, str] = {}
        self._prepare_maps()

    def _load_catalog(self) -> List[Dict[str, Any]]:
        with open(self.catalog_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)["modules"]

    def _prepare_maps(self):
        for m in self.modules_cfg:
            m_id = m["id"]
            for path_str in m["paths"]:
                p = (self.repo_root / path_str).resolve()
                if p.is_dir():
                    for py_file in p.rglob("*.py"):
                        if "__pycache__" in py_file.parts or py_file.suffix == ".pyc":
                            continue
                        fp = py_file.resolve()
                        if fp in self.file_to_module:
                            raise ValueError(f"File {fp} belongs to multiple modules")
                        self.file_to_module[fp] = m_id
                else:
                    fp = p.resolve()
                    if fp in self.file_to_module:
                        raise ValueError(f"File {fp} belongs to multiple modules")
                    self.file_to_module[fp] = m_id

        src_chemuson = self.repo_root / "src" / "chemuson"
        for fp, m_id in self.file_to_module.items():
            if src_chemuson in fp.parents or fp == src_chemuson:
                try:
                    rel_path = fp.relative_to(src_chemuson)
                    parts = list(rel_path.parts)
                    if parts[-1] == "__init__.py":
                        parts.pop()
                    elif parts[-1].endswith(".py"):
                        parts[-1] = parts[-1][:-3]

                    if not parts or parts == ["__init__"]:
                        mod_name = "chemuson"
                    else:
                        mod_name = "chemuson." + ".".join(parts)
                    self.importable_to_module_id[mod_name] = m_id
                except ValueError:
                    continue

    def _get_module_name(self, fp: Path) -> str:
        src_chemuson = self.repo_root / "src" / "chemuson"
        if not (fp.resolve().is_relative_to(src_chemuson.resolve())):
            raise ValueError(
                f"Path {fp} is not under {src_chemuson}"
            )
        rel_path = fp.relative_to(src_chemuson)
        parts = list(rel_path.parts)
        if parts[-1] == "__init__.py":
            parts.pop()
        elif parts[-1].endswith(".py"):
            parts[-1] = parts[-1][:-3]
        if not parts or parts == ["__init__"]:
            return "chemuson"
        return "chemuson." + ".".join(parts)

    def _resolve_import(self, imp_data: Dict[str, Any], fp: Path) -> tuple[str, bool] | None:
        target_mod_name = ""
        if imp_data["type"] == "absolute":
            target_mod_name = imp_data["name"]
        elif imp_data["type"] == "from":
            level = imp_data["level"]
            if level == 0:
                target_mod_name = imp_data["module"]
            else:
                current_mod = self._get_module_name(fp)
                parts = current_mod.split(".")
                if fp.name == "__init__.py":
                    base_parts = parts[:len(parts) - level + 1]
                else:
                    base_parts = parts[:-level]

                if not base_parts:
                    base_parts = ["chemuson"]

                if imp_data["module"]:
                    target_mod_name = ".".join(base_parts + [imp_data["module"]])
                else:
                    target_mod_name = ".".join(base_parts)

        if not target_mod_name:
            return None

        best_match_id = None
        best_match_len = -1
        for imp_mod, m_id in self.importable_to_module_id.items():
            if target_mod_name == imp_mod or target_mod_name.startswith(imp_mod + "."):
                if len(imp_mod) > best_match_len:
                    best_match_len = len(imp_mod)
                    best_match_id = m_id

        if best_match_id:
            return (best_match_id, imp_data["type_checking_only"])
        return None

    def analyze_all(self) -> List[ObservedImport]:
        all_observed = []
        for fp, m_id in self.file_to_module.items():
            if not fp.exists():
                continue
            with open(fp, "r", encoding="utf-8") as f:
                source = f.read()
            try:
                tree = ast.parse(source)
            except SyntaxError as exc:
                raise AssertionError(
                    f"Cannot parse {fp} at line {exc.lineno}: {exc.msg}"
                ) from exc

            visitor = ImportVisitor(source, fp, self._get_module_name(fp))
            visitor.visit(tree)

            for imp_data in visitor.imports:
                resolved = self._resolve_import(imp_data, fp)
                if resolved:
                    target_id, tc_only = resolved
                    if target_id != m_id:
                        node = imp_data["node"]
                        stmt = ast.get_source_segment(source, node)
                        if stmt is None:
                            stmt = ast.unparse(node)
                        stmt = re.sub(r"\s+", " ", stmt).strip()

                        all_observed.append(ObservedImport(
                            source_id=m_id,
                            target_id=target_id,
                            file=fp,
                            line=imp_data["lineno"],
                            statement=stmt,
                            type_checking_only=tc_only
                        ))
        return all_observed


# ---------------------------------------------------------------------------
# Helper functions for path and import normalization
# ---------------------------------------------------------------------------

def _normalize_policy_path(repo_root: Path, value: str | Path) -> str:
    """
    Normalize a path for policy matching.
    - Accepts str or Path.
    - Converts all backslashes to forward slashes first.
    - If the path is absolute and under repo_root, returns the relative POSIX path.
    - If the path is absolute but outside repo_root, returns the normalized absolute POSIX path.
    - If the path is relative, returns the normalized relative POSIX path.
    - Does not require the file to exist.
    - Does not use resolve() for synthetic files.
    - Does not depend on current working directory.
    """
    # First, normalize backslashes to forward slashes in the string representation
    # We need to do this carefully because Path may interpret \ as part of the name on POSIX
    # So we replace in the string before creating Path, but after ensuring it's a string
    if isinstance(value, str):
        value_str = value.replace("\\", "/")
    else:
        value_str = str(value).replace("\\", "/")

    # Recreate Path from the normalized string
    p = Path(value_str)

    # Check if it's an absolute path under repo_root
    try:
        rel = p.relative_to(repo_root)
        return str(rel).replace("\\", "/")
    except ValueError:
        # Not under repo_root, return as normalized absolute or relative path
        # Use as_posix() to ensure forward slashes
        return p.as_posix()


def _normalize_import_path(value: str) -> str:
    """Normalize import statement by collapsing whitespace and stripping."""
    return re.sub(r"\s+", " ", value).strip()


# ---------------------------------------------------------------------------
# TYPE_CHECKING Exception Validation Helper
# ---------------------------------------------------------------------------

def check_type_checking_exceptions(
    analyzer: ImportBoundaryAnalyzer,
    catalog: list,
    observed_imports: list,
) -> list[dict]:
    """
    Check documented TYPE_CHECKING exceptions against observed TYPE_CHECKING imports.
    Returns a list of exceptions that are documented but do not correspond exactly
    to any observed TYPE_CHECKING import.

    Identity uses: source_id, target_id, normalized file, normalized import_path, True.
    Only imports with type_checking_only=True can satisfy a TYPE_CHECKING exception.
    An observed TYPE_CHECKING import without a documented exception is allowed.
    """
    # Build set of documented TYPE_CHECKING exception keys and mapping
    exception_keys = set()
    exceptions_by_key = {}

    for m in catalog:
        for exc in m.get("temporary_exceptions", []):
            if exc.get("type_checking_only") is True:
                source_id = exc["source_id"]
                file_posix = _normalize_policy_path(analyzer.repo_root, exc["file"])
                import_path_norm = _normalize_import_path(exc["import_path"])
                key = (source_id, exc["target_id"], file_posix, import_path_norm, True)
                exception_keys.add(key)
                exceptions_by_key[key] = exc

    # Build set of observed TYPE_CHECKING import keys
    observed_keys = set()
    for imp in observed_imports:
        if imp.type_checking_only:
            source_id = imp.source_id
            target_id = imp.target_id
            file_posix = _normalize_policy_path(analyzer.repo_root, imp.file)
            import_path_norm = _normalize_import_path(imp.statement)
            key = (source_id, target_id, file_posix, import_path_norm, True)
            observed_keys.add(key)

    # Find documented exceptions that have no matching observed import
    obsolete_exceptions = []
    for key in exception_keys:
        if key not in observed_keys:
            exc = exceptions_by_key[key]
            obsolete_exceptions.append(exc)

    return obsolete_exceptions


# ---------------------------------------------------------------------------
# Policy Enforcement Helper
# ---------------------------------------------------------------------------

def check_runtime_policy(analyzer: ImportBoundaryAnalyzer, catalog: list, observed_imports: list):
    """
    Check runtime import policy for observed imports.
    Returns (violations, obsolete_exceptions) where each is a list of dicts.
    Violations have keys: file, line, source_id, target_id, statement, is_forbidden, not_in_target.
    Obsolete exceptions have the exception data from the catalog.
    """
    modules_by_id = {m["id"]: m for m in catalog}

    # Build exception key set and mapping for runtime (type_checking_only=False) exceptions
    exception_keys = set()
    exceptions_by_key = {}

    for m in catalog:
        for exc in m.get("temporary_exceptions", []):
            if exc.get("type_checking_only") is False:
                # Use the source_id from the exception itself, not the module id
                source_id = exc["source_id"]
                # Normalize file to POSIX relative path using helper
                file_posix = _normalize_policy_path(analyzer.repo_root, exc["file"])
                import_path_norm = _normalize_import_path(exc["import_path"])
                key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                exception_keys.add(key)
                exceptions_by_key[key] = exc

    violations = []
    matched_exceptions = set()

    for imp in observed_imports:
        if imp.type_checking_only:
            continue  # Skip TYPE_CHECKING imports for runtime policy

        source_id = imp.source_id
        target_id = imp.target_id
        m_config = modules_by_id.get(source_id)
        if not m_config:
            continue

        # Normalize observed file to POSIX relative path using helper
        file_posix = _normalize_policy_path(analyzer.repo_root, imp.file)
        import_path_norm = _normalize_import_path(imp.statement)

        target_deps = set(m_config.get("target_dependencies", []))
        forbidden_deps = set(m_config.get("forbidden_dependencies", []))

        is_forbidden = target_id in forbidden_deps
        not_in_target = target_id not in target_deps

        if is_forbidden or not_in_target:
            exc_key = (source_id, target_id, file_posix, import_path_norm, False)
            if exc_key in exception_keys:
                matched_exceptions.add(exc_key)
                continue
            # Violation
            violations.append({
                "file": imp.file,
                "line": imp.line,
                "source_id": source_id,
                "target_id": target_id,
                "statement": imp.statement,
                "is_forbidden": is_forbidden,
                "not_in_target": not_in_target,
            })

    # Check for obsolete exceptions
    obsolete_exceptions = []
    for key in exception_keys:
        if key not in matched_exceptions:
            exc = exceptions_by_key[key]
            obsolete_exceptions.append(exc)

    return violations, obsolete_exceptions


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.fixture
def analyzer():
    repo_root = Path(__file__).resolve().parent.parent.parent
    catalog_path = repo_root / "architecture" / "modules.yml"
    return ImportBoundaryAnalyzer(repo_root, catalog_path)


class TestImportBoundaries:

    def test_module_file_ownership_has_no_overlap(self, analyzer):
        assert len(analyzer.file_to_module) > 0

    def test_relative_import_resolution_examples(self, analyzer):
        src_chemuson = analyzer.repo_root / "src" / "chemuson"
        expected_ids = {module["id"] for module in analyzer.modules_cfg}

        cases = [
            {
                "path": src_chemuson / "chemcalc" / "formula.py",
                "level": 1,
                "module": "valence",
                "expected_target": "M03",
                "description": "chemcalc/formula.py: from .valence import implicit_h_count"
            },
            {
                "path": src_chemuson / "chemname" / "coordination.py",
                "level": 1,
                "module": "molview",
                "expected_target": "M04",
                "description": "chemname/coordination.py: from .molview import MolView"
            },
            {
                "path": src_chemuson / "gui" / "canvas" / "canvas_view.py",
                "level": 1,
                "module": "canvas_constants",
                "expected_target": "M09",
                "description": "gui/canvas/canvas_view.py: from .canvas_constants import ..."
            },
        ]

        for case in cases:
            fp = case["path"]
            with open(fp, "r", encoding="utf-8") as f:
                source_text = f.read()
            tree = ast.parse(source_text)

            imp_node = None
            for node in ast.walk(tree):
                if (isinstance(node, ast.ImportFrom) and
                        node.level == case["level"] and
                        node.module == case["module"]):
                    imp_node = node
                    break

            assert imp_node is not None, (
                f"Expected import not found in {fp}: "
                f"level={case['level']}, module={case['module']} "
                f"({case['description']})"
            )

            m_id = analyzer.file_to_module[fp]
            assert m_id in expected_ids, f"source_id '{m_id}' not in catalog"

            imp_data = {
                "type": "from",
                "module": imp_node.module,
                "level": imp_node.level,
                "names": [a.name for a in imp_node.names],
                "lineno": imp_node.lineno,
                "node": imp_node,
                "type_checking_only": False
            }
            resolved = analyzer._resolve_import(imp_data, fp)
            assert resolved is not None, (
                f"Could not resolve relative import in {fp} at line {imp_node.lineno}"
            )
            target_id, tc_only = resolved
            assert target_id == case["expected_target"], (
                f"Expected target {case['expected_target']} but got {target_id} "
                f"for {case['description']}"
            )
            assert tc_only is False, (
                f"Expected type_checking_only=False for {case['description']}"
            )

    def test_type_checking_imports_are_classified_separately(self, analyzer):
        src_chemuson = analyzer.repo_root / "src" / "chemuson"
        persistence_file = src_chemuson / "chemio" / "persistence.py"
        autosave_file = src_chemuson / "utils" / "autosave.py"

        all_imports = analyzer.analyze_all()

        p_imports = [i for i in all_imports if i.file == persistence_file]
        canvas_tc = [i for i in p_imports if i.target_id == "M09" and "ChemusonCanvas" in i.statement]
        assert len(canvas_tc) > 0
        assert canvas_tc[0].type_checking_only is True

        a_imports = [i for i in all_imports if i.file == autosave_file]
        canvas_tc_a = [i for i in a_imports if i.target_id == "M09" and "ChemusonCanvas" in i.statement]
        assert len(canvas_tc_a) > 0
        assert canvas_tc_a[0].type_checking_only is True

    def test_type_checking_else_block_is_runtime(self, analyzer):
        synthetic_source = (
            "if TYPE_CHECKING:\n"
            "    from chemuson.gui.canvas import ChemusonCanvas\n"
            "else:\n"
            "    from chemuson.core import MolGraph\n"
        )
        synthetic_tree = ast.parse(synthetic_source)
        synthetic_visitor = ImportVisitor(
            synthetic_source, Path("/tmp/fake.py"), "chemuson.fake"
        )
        synthetic_visitor.visit(synthetic_tree)

        tc_imports = [imp for imp in synthetic_visitor.imports if imp["type_checking_only"]]
        rt_imports = [imp for imp in synthetic_visitor.imports if not imp["type_checking_only"]]

        assert len(tc_imports) == 1, (
            f"Expected exactly 1 TYPE_CHECKING import, found {len(tc_imports)}"
        )
        assert "ChemusonCanvas" in tc_imports[0]["names"][0], (
            f"Expected ChemusonCanvas in TYPE_CHECKING import, got {tc_imports[0]['names']}"
        )

        assert len(rt_imports) == 1, (
            f"Expected exactly 1 runtime import in else block, found {len(rt_imports)}"
        )
        assert "MolGraph" in rt_imports[0]["names"][0], (
            f"Expected MolGraph in runtime import, got {rt_imports[0]['names']}"
        )

        # Resolver los imports a través del analyzer para verificar mapeo y booleanos
        imp_data_tc = tc_imports[0]
        imp_data_rt = rt_imports[0]

        # Para el import de TYPE_CHECKING, debe resolverse a M09 (gui.canvas)
        resolved_tc = analyzer._resolve_import(imp_data_tc, Path("/tmp/fake.py"))
        assert resolved_tc is not None, (
            f"Failed to resolve TYPE_CHECKING import {imp_data_tc}"
        )
        target_id_tc, tc_only_flag = resolved_tc
        assert target_id_tc == "M09", (
            f"Expected M09 for TYPE_CHECKING import, got {target_id_tc}"
        )
        assert tc_only_flag is True, (
            f"Expected type_checking_only=True for TYPE_CHECKING import, got {tc_only_flag}"
        )

        # Para el import del else, debe resolverse a M00 (core)
        resolved_rt = analyzer._resolve_import(imp_data_rt, Path("/tmp/fake.py"))
        assert resolved_rt is not None, (
            f"Failed to resolve runtime import {imp_data_rt}"
        )
        target_id_rt, tc_only_flag_rt = resolved_rt
        assert target_id_rt == "M00", (
            f"Expected M00 for runtime import, got {target_id_rt}"
        )
        assert tc_only_flag_rt is False, (
            f"Expected type_checking_only=False for runtime import, got {tc_only_flag_rt}"
        )

    def test_typing_type_checking_import_classification(self, analyzer):
        """Sintético para verificar typing.TYPE_CHECKING y else."""
        synthetic_source = (
            "if typing.TYPE_CHECKING:\n"
            "    from chemuson.gui.canvas import ChemusonCanvas\n"
            "else:\n"
            "    from chemuson.core import MolGraph\n"
        )
        synthetic_tree = ast.parse(synthetic_source)
        synthetic_visitor = ImportVisitor(
            synthetic_source, Path("/tmp/fake.py"), "chemuson.fake"
        )
        synthetic_visitor.visit(synthetic_tree)

        tc_imports = [imp for imp in synthetic_visitor.imports if imp["type_checking_only"]]
        rt_imports = [imp for imp in synthetic_visitor.imports if not imp["type_checking_only"]]

        assert len(tc_imports) == 1, (
            f"Expected exactly 1 TYPE_CHECKING import, found {len(tc_imports)}"
        )
        assert "ChemusonCanvas" in tc_imports[0]["names"][0], (
            f"Expected ChemusonCanvas in TYPE_CHECKING import, got {tc_imports[0]['names']}"
        )

        assert len(rt_imports) == 1, (
            f"Expected exactly 1 runtime import in else block, found {len(rt_imports)}"
        )
        assert "MolGraph" in rt_imports[0]["names"][0], (
            f"Expected MolGraph in runtime import, got {rt_imports[0]['names']}"
        )

        imp_data_tc = tc_imports[0]
        imp_data_rt = rt_imports[0]

        resolved_tc = analyzer._resolve_import(imp_data_tc, Path("/tmp/fake.py"))
        assert resolved_tc is not None, (
            f"Failed to resolve TYPE_CHECKING import {imp_data_tc}"
        )
        target_id_tc, tc_only_flag = resolved_tc
        assert target_id_tc == "M09", (
            f"Expected M09 for typing.TYPE_CHECKING import, got {target_id_tc}"
        )
        assert tc_only_flag is True, (
            f"Expected type_checking_only=True for typing.TYPE_CHECKING import, got {tc_only_flag}"
        )

        resolved_rt = analyzer._resolve_import(imp_data_rt, Path("/tmp/fake.py"))
        assert resolved_rt is not None, (
            f"Failed to resolve runtime import {imp_data_rt}"
        )
        target_id_rt, tc_only_flag_rt = resolved_rt
        assert target_id_rt == "M00", (
            f"Expected M00 for runtime import, got {target_id_rt}"
        )
        assert tc_only_flag_rt is False, (
            f"Expected type_checking_only=False for runtime import, got {tc_only_flag_rt}"
        )

    def test_runtime_dependencies_match_catalog(self, analyzer):
        all_observed = analyzer.analyze_all()

        for m in analyzer.modules_cfg:
            m_id = m["id"]
            expected = set(m["current_dependencies"])

            observed_for_m = set()
            for imp in all_observed:
                if not imp.type_checking_only and analyzer.file_to_module.get(imp.file) == m_id:
                    observed_for_m.add(imp.target_id)

            expected = {e for e in expected if e != m_id}

            missing_in_catalog = observed_for_m - expected
            declared_but_not_observed = expected - observed_for_m

            assert missing_in_catalog == set(), (
                f"Module {m_id} ({m['name']}) has runtime dependencies not in catalog: {missing_in_catalog}"
            )
            assert declared_but_not_observed == set(), (
                f"Module {m_id} ({m['name']}) has declared dependencies not observed at runtime: {declared_but_not_observed}"
            )

    def test_observed_imports_have_valid_metadata(self, analyzer):
        all_observed = analyzer.analyze_all()
        expected_ids = {module["id"] for module in analyzer.modules_cfg}
        for imp in all_observed:
            assert imp.source_id in expected_ids, (
                f"source_id '{imp.source_id}' not in catalog IDs"
            )
            assert imp.target_id in expected_ids, (
                f"target_id '{imp.target_id}' not in catalog IDs"
            )
            assert imp.source_id != imp.target_id, (
                f"source_id equals target_id in {imp.file}:{imp.line}"
            )
            assert imp.file.exists(), f"File does not exist: {imp.file}"
            assert imp.line > 0, f"Line must be positive: {imp.line}"
            assert len(imp.statement) > 0, f"Statement is empty in {imp.file}:{imp.line}"
            assert re.match(r"^M\d\d$", imp.target_id), (
                f"target_id '{imp.target_id}' does not match ^M\\d\\d$"
            )

    def test_runtime_policy_enforcement(self, analyzer):
        """
        Enforce runtime import policy on the actual catalog and codebase.
        This test verifies that all observed runtime imports either:
        - are allowed by target_dependencies/forbidden_dependencies; or
        - have an exact temporary_exception matching source_id, target_id, file, import_path, type_checking_only.
        Also verifies there are no obsolete runtime exceptions.
        """
        all_observed = analyzer.analyze_all()
        violations, obsolete_exceptions = check_runtime_policy(analyzer, analyzer.modules_cfg, all_observed)

        for v in violations:
            rules = []
            if v["is_forbidden"]:
                rules.append("forbidden_dependencies")
            if v["not_in_target"]:
                rules.append("target_dependencies")
            rule_str = " and ".join(rules)
            assert False, (
                f"Runtime boundary violation in {v['file']}:{v['line']} | "
                f"Module {v['source_id']} -> {v['target_id']} | "
                f"Import: {v['statement']} | "
                f"Violates: {rule_str} | "
                f"No exact temporary_exception found"
            )

        for exc in obsolete_exceptions:
            assert False, (
                f"Obsolete runtime exception: source={exc['source_id']}, target={exc['target_id']}, "
                f"file={exc['file']}, import_path={exc['import_path']}, type_checking_only=False | "
                f"This exception does not correspond to any observed runtime violation"
            )

    def test_type_checking_exceptions_documented_are_valid(self, analyzer):
        """
        Validate that all TYPE_CHECKING exceptions in the real catalog correspond
        to actual observed TYPE_CHECKING imports. No obsolete TYPE_CHECKING exceptions.
        """
        all_observed = analyzer.analyze_all()
        obsolete = check_type_checking_exceptions(
            analyzer,
            analyzer.modules_cfg,
            all_observed,
        )
        if obsolete:
            msg_lines = [f"Found {len(obsolete)} obsolete TYPE_CHECKING exceptions. Each documented TYPE_CHECKING exception must match an observed import."]
            for exc in obsolete:
                msg_lines.append(
                    f"  - source_id={exc['source_id']}, target_id={exc['target_id']}, "
                    f"file={exc['file']}, import_path={exc['import_path']}, "
                    f"type_checking_only={exc['type_checking_only']} | "
                    f"Este excepción no corresponde a ningún import TYPE_CHECKING observado."
                )
            assert False, "\n".join(msg_lines)


    # Synthetic tests that directly call check_runtime_policy

    def test_synthetic_allowed_dependency(self, analyzer):
        """Case 1: Dependency in target_dependencies and not forbidden -> no violations."""
        test_catalog = [
            {"id": "M99", "target_dependencies": ["M00"], "forbidden_dependencies": [], "temporary_exceptions": []},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M00", file=Path("/tmp/test.py"), line=1,
                           statement="from chemuson.core import MolGraph", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_missing_target_dependency(self, analyzer):
        """Case 2: Dependency not in target_dependencies -> violation."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M01", file=Path("/tmp/test.py"), line=2,
                           statement="from chemuson.chemio import something", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert violations[0]["not_in_target"] is True
        assert violations[0]["is_forbidden"] is False

    def test_synthetic_forbidden_dependency(self, analyzer):
        """Case 3: Dependency in forbidden_dependencies -> violation."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"], "temporary_exceptions": []},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert violations[0]["is_forbidden"] is True

    def test_synthetic_both_forbidden_and_missing(self, analyzer):
        """Case 4: Dependency both forbidden and missing from target -> violation with both flags."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"], "temporary_exceptions": []},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert violations[0]["is_forbidden"] is True
        assert violations[0]["not_in_target"] is True

    def test_synthetic_exact_exception_allows(self, analyzer):
        """Case 5: Exact exception allows violation -> no violation, exception not obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_wrong_source_id_exception(self, analyzer):
        """Case 6: Exception with wrong source_id -> violation, exception obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M98", "target_id": "M02", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert len(obsolete) == 1
        assert obsolete[0]["source_id"] == "M98"

    def test_synthetic_wrong_target_id_exception(self, analyzer):
        """Case 7: Exception with wrong target_id -> violation, exception obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M03", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert len(obsolete) == 1

    def test_synthetic_wrong_file_exception(self, analyzer):
        """Case 8: Exception with wrong file -> violation, exception obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02", "file": "/tmp/wrong.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert len(obsolete) == 1

    def test_synthetic_wrong_import_path_exception(self, analyzer):
        """Case 9: Exception with wrong import_path -> violation, exception obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import wrong",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert len(obsolete) == 1

    def test_synthetic_type_checking_exception_on_runtime_import(self, analyzer):
        """Case 10: type_checking_only=True exception does not allow runtime import."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 1
        assert len(obsolete) == 0

    def test_synthetic_exception_without_observed_import(self, analyzer):
        """Case 11: Runtime exception with no observed import -> obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": ["M00"], "forbidden_dependencies": [],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M00", file=Path("/tmp/test.py"), line=1,
                           statement="from chemuson.core import MolGraph", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 1

    def test_synthetic_exception_on_allowed_import_is_obsolete(self, analyzer):
        """Case 12: Exception that corresponds to an already allowed import -> obsolete."""
        test_catalog = [
            {"id": "M99", "target_dependencies": ["M00"], "forbidden_dependencies": [],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M00", "file": "/tmp/test.py",
                                      "import_path": "from chemuson.core import MolGraph",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M00", file=Path("/tmp/test.py"), line=1,
                           statement="from chemuson.core import MolGraph", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 1

    def test_synthetic_relative_vs_absolute_path(self, analyzer):
        """Case 13a: Relative path from catalog matches absolute path observed under repo_root."""
        repo_file = Path("tests/architecture/test_import_boundaries.py").resolve()
        assert analyzer.repo_root in repo_file.parents or analyzer.repo_root == repo_file
        rel_path = str(repo_file.relative_to(analyzer.repo_root))
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": rel_path,
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=repo_file, line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_backslash_vs_forward_slash(self, analyzer):
        """Case 13b: Backslash in path string should normalize to forward slash."""
        # Simulate a path with backslashes (as might come from Windows or raw strings)
        # The file itself is under repo_root
        repo_file = Path("tests/architecture/test_import_boundaries.py").resolve()
        assert analyzer.repo_root in repo_file.parents or analyzer.repo_root == repo_file

        # Construct a string with backslashes
        backslash_path = str(repo_file.relative_to(analyzer.repo_root)).replace("/", "\\")
        # This should be something like "tests\\architecture\\test_import_boundaries.py"

        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": backslash_path,
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=repo_file, line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_whitespace_equivalence(self, analyzer):
        """Case 13c: Whitespace variations in import_path should match."""
        repo_file = Path("tests/architecture/test_import_boundaries.py").resolve()
        assert analyzer.repo_root in repo_file.parents or analyzer.repo_root == repo_file
        rel_path = str(repo_file.relative_to(analyzer.repo_root))

        # Import statement with extra whitespace, tabs, newlines
        import_with_extra_ws = "from  chemuson.clean2d  import\tlayout\n"

        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": rel_path,
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=repo_file, line=3,
                           statement=import_with_extra_ws, type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_external_path_stability(self, analyzer):
        """Case 13d: External paths should normalize stably without becoming relative to repo."""
        # Use an external absolute path (not under repo_root)
        external_path = Path("/some/external/path/file.py")

        # Normalized form should be the same regardless of how it's represented
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": [],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": str(external_path),
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        # Observe the same external path but as a Path object
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=external_path, line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        # Should match because both normalize to same POSIX string
        assert len(violations) == 0
        assert len(obsolete) == 0

    def test_synthetic_real_difference_fails(self, analyzer):
        """Case 13e: Real path differences should still fail to match."""
        repo_file = Path("tests/architecture/test_import_boundaries.py").resolve()
        assert analyzer.repo_root in repo_file.parents or analyzer.repo_root == repo_file
        rel_path = str(repo_file.relative_to(analyzer.repo_root))

        # Use a different file path
        different_path = "tests/architecture/test_module_catalog.py"

        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": different_path,
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=repo_file, line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        violations, obsolete = check_runtime_policy(analyzer, test_catalog, observed)
        # Should fail to match because paths are different
        assert len(violations) == 1
        assert len(obsolete) == 1

    def test_synthetic_type_checking_exception_exact_match(self, analyzer):
        """Case 1: Excepción TYPE_CHECKING exacta -> no obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 0

    def test_synthetic_type_checking_wrong_source_id(self, analyzer):
        """Case 2: source_id incorrecto -> obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M98", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 1
        assert obsolete[0]["source_id"] == "M98"

    def test_synthetic_type_checking_wrong_target_id(self, analyzer):
        """Case 3: target_id incorrecto -> obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M03",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 1

    def test_synthetic_type_checking_wrong_file(self, analyzer):
        """Case 4: archivo incorrecto -> obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/wrong.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 1

    def test_synthetic_type_checking_wrong_import_path(self, analyzer):
        """Case 5: import_path incorrecto -> obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import wrong",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 1

    def test_synthetic_type_checking_exception_vs_runtime_import(self, analyzer):
        """Case 6: Excepción TYPE_CHECKING frente a import runtime -> obsoleta."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 1

    def test_synthetic_type_checking_runtime_exception_ignored(self, analyzer):
        """Case 7: Excepción runtime frente a import TYPE_CHECKING -> se ignora en esta función."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": False, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 0

    def test_synthetic_type_checking_import_without_exception_allowed(self, analyzer):
        """Case 8: Import TYPE_CHECKING sin excepción documentada -> no produce error."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": []},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 0

    def test_synthetic_type_checking_normalization(self, analyzer):
        """Case 9: Normalización - ruta relativa vs absoluta, backslashes, whitespace."""
        repo_file = Path("tests/architecture/test_import_boundaries.py").resolve()
        assert analyzer.repo_root in repo_file.parents or analyzer.repo_root == repo_file
        rel_path = str(repo_file.relative_to(analyzer.repo_root))

        # Import with extra whitespace and backslashes in path string
        import_with_extra_ws = "from  chemuson.clean2d  import\tlayout\n"
        backslash_path = rel_path.replace("/", "\\")

        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": backslash_path,
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=repo_file, line=3,
                           statement=import_with_extra_ws, type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 0

    def test_synthetic_type_checking_two_imports_one_tc_one_rt(self, analyzer):
        """Case 10: Dos imports observados, uno runtime y otro TYPE_CHECKING -> solo TC satisface."""
        test_catalog = [
            {"id": "M99", "target_dependencies": [], "forbidden_dependencies": ["M02"],
             "temporary_exceptions": [{"source_id": "M99", "target_id": "M02",
                                      "file": "/tmp/test.py",
                                      "import_path": "from chemuson.clean2d import layout",
                                      "type_checking_only": True, "reason": "test", "debt_ref": "test",
                                      "elimination_condition": "fix it"}]},
            {"id": "M00", "target_dependencies": [], "forbidden_dependencies": [], "temporary_exceptions": []},
        ]
        observed = [
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=False),
            ObservedImport(source_id="M99", target_id="M02", file=Path("/tmp/test.py"), line=3,
                           statement="from chemuson.clean2d import layout", type_checking_only=True)
        ]
        obsolete = check_type_checking_exceptions(analyzer, test_catalog, observed)
        assert len(obsolete) == 0
