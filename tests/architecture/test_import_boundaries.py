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
        source_id = m["id"]
        for exc in m.get("temporary_exceptions", []):
            if exc.get("type_checking_only") is False:
                # Normalize file to POSIX relative path
                try:
                    file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    # If file not under repo root, use absolute as is but normalized
                    file_posix = str(Path(exc["file"]).as_posix())
                import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
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

        # Normalize observed file to POSIX relative path
        try:
            file_rel = imp.file.relative_to(analyzer.repo_root)
            file_posix = str(file_rel).replace("\\", "/")
        except ValueError:
            file_posix = str(imp.file.as_posix())

        import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()

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

    def test_policy_enforcement_synthetic_cases(self, analyzer):
        """
        Synthetic tests for policy enforcement:
        1. Allowed dependency (in target_dependencies, not forbidden)
        2. Dependency not in target_dependencies
        3. Dependency in forbidden_dependencies
        4. Exception exact match allows violation
        5. Exception with wrong file does NOT allow
        6. Exception with wrong import_path does NOT allow
        7. type_checking_only=True exception does NOT allow runtime import
        8. Runtime exception without observed violation is obsolete
        """
        modules_by_id = {m["id"]: m for m in analyzer.modules_cfg}

        # Build a temporary catalog with controlled rules for testing
        test_catalog = [
            {
                "id": "M99",
                "name": "test_module",
                "target_dependencies": ["M00"],
                "forbidden_dependencies": [],
                "temporary_exceptions": []
            },
            {
                "id": "M00",
                "name": "core",
                "target_dependencies": [],
                "forbidden_dependencies": [],
                "temporary_exceptions": []
            }
        ]

        # Test 1: Allowed dependency (no violation)
        observed_allowed = [
            ObservedImport(
                source_id="M99",
                target_id="M00",
                file=Path("/tmp/test.py"),
                line=1,
                statement="from chemuson.core import MolGraph",
                type_checking_only=False
            )
        ]
        # This should pass without assertion because M00 is in target_dependencies

        # We'll actually call enforce_runtime_policy logic here but with synthetic data.
        # Instead of duplicating, we test by manipulating analyzer state minimally.
        # Better: create a separate function to call. But the spec says don't introduce new architecture.
        # So we inline the check.

        # Build exception keys for this synthetic catalog
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)

        # Simulate enforcement for allowed case
        violations = []
        for imp in observed_allowed:
            m_config = modules_by_id.get(imp.source_id) or next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                # Normalize file path similarly to check_runtime_policy
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 0, "Allowed import should not produce violation"

        # Test 2: Dependency not in target_dependencies -> violation
        observed_not_in_target = [
            ObservedImport(
                source_id="M99",
                target_id="M01",
                file=Path("/tmp/test.py"),
                line=2,
                statement="from chemuson.chemio import something",
                type_checking_only=False
            )
        ]
        violations = []
        for imp in observed_not_in_target:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 1, "Should report violation for missing target dependency"
        assert violations[0]["not_in_target"] is True
        assert violations[0]["is_forbidden"] is False

        # Test 3: Dependency in forbidden_dependencies -> violation
        # Update test_catalog to forbid M02
        test_catalog[0]["forbidden_dependencies"] = ["M02"]
        observed_forbidden = [
            ObservedImport(
                source_id="M99",
                target_id="M02",
                file=Path("/tmp/test.py"),
                line=3,
                statement="from chemuson.clean2d import layout",
                type_checking_only=False
            )
        ]
        violations = []
        for imp in observed_forbidden:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 1, "Should report violation for forbidden dependency"
        assert violations[0]["is_forbidden"] is True

        # Test 4: Exact exception allows violation
        test_catalog[0]["temporary_exceptions"] = [
            {
                "source_id": "M99",
                "target_id": "M02",
                "file": "/tmp/test.py",
                "import_path": "from chemuson.clean2d import layout",
                "type_checking_only": False,
                "reason": "test exception",
                "debt_ref": "test",
                "elimination_condition": "fix it"
            }
        ]
        # Rebuild exception keys
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)
        observed_with_exception = [
            ObservedImport(
                source_id="M99",
                target_id="M02",
                file=Path("/tmp/test.py"),
                line=3,
                statement="from chemuson.clean2d import layout",
                type_checking_only=False
            )
        ]
        violations = []
        for imp in observed_with_exception:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 0, "Exact exception should allow violation"

        # Test 5: Exception with wrong file does NOT allow
        test_catalog[0]["temporary_exceptions"] = [
            {
                "source_id": "M99",
                "target_id": "M02",
                "file": "/tmp/wrong_file.py",
                "import_path": "from chemuson.clean2d import layout",
                "type_checking_only": False,
                "reason": "test exception",
                "debt_ref": "test",
                "elimination_condition": "fix it"
            }
        ]
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)
        violations = []
        for imp in observed_with_exception:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 1, "Wrong file should not allow violation"

        # Test 6: Exception with wrong import_path does NOT allow
        test_catalog[0]["temporary_exceptions"] = [
            {
                "source_id": "M99",
                "target_id": "M02",
                "file": "/tmp/test.py",
                "import_path": "from chemuson.clean2d import something_else",
                "type_checking_only": False,
                "reason": "test exception",
                "debt_ref": "test",
                "elimination_condition": "fix it"
            }
        ]
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)
        violations = []
        for imp in observed_with_exception:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 1, "Wrong import_path should not allow violation"

        # Test 7: type_checking_only=True exception does NOT allow runtime import
        test_catalog[0]["temporary_exceptions"] = [
            {
                "source_id": "M99",
                "target_id": "M02",
                "file": "/tmp/test.py",
                "import_path": "from chemuson.clean2d import layout",
                "type_checking_only": True,
                "reason": "test exception",
                "debt_ref": "test",
                "elimination_condition": "fix it"
            }
        ]
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)
        # Note: we are NOT adding the True exception to exception_keys because it's type_checking_only=True
        violations = []
        for imp in observed_with_exception:
            m_config = next((m for m in test_catalog if m["id"] == imp.source_id), None)
            if not m_config:
                continue
            target_deps = set(m_config.get("target_dependencies", []))
            forbidden_deps = set(m_config.get("forbidden_dependencies", []))
            is_forbidden = imp.target_id in forbidden_deps
            not_in_target = imp.target_id not in target_deps
            if is_forbidden or not_in_target:
                try:
                    file_rel = imp.file.relative_to(analyzer.repo_root)
                    file_posix = str(file_rel).replace("\\", "/")
                except ValueError:
                    file_posix = str(imp.file.as_posix())
                import_path_norm = re.sub(r"\s+", " ", imp.statement).strip()
                exc_key = (imp.source_id, imp.target_id, file_posix, import_path_norm, False)
                if exc_key not in exception_keys:
                    violations.append({
                        "file": imp.file,
                        "line": imp.line,
                        "source_id": imp.source_id,
                        "target_id": imp.target_id,
                        "statement": imp.statement,
                        "is_forbidden": is_forbidden,
                        "not_in_target": not_in_target,
                    })
        assert len(violations) == 1, "type_checking_only=True exception should not allow runtime import"

        # Test 8: Runtime exception without observed violation is obsolete
        test_catalog[0]["temporary_exceptions"] = [
            {
                "source_id": "M99",
                "target_id": "M02",
                "file": "/tmp/test.py",
                "import_path": "from chemuson.clean2d import layout",
                "type_checking_only": False,
                "reason": "test exception",
                "debt_ref": "test",
                "elimination_condition": "fix it"
            }
        ]
        # No observed import that matches this exception
        # We need to simulate the obsolete check
        matched_exceptions = set()  # Empty because no matching observed import
        exception_keys = set()
        for m in test_catalog:
            source_id = m["id"]
            for exc in m.get("temporary_exceptions", []):
                if exc.get("type_checking_only") is False:
                    try:
                        file_rel = Path(exc["file"]).relative_to(analyzer.repo_root)
                        file_posix = str(file_rel).replace("\\", "/")
                    except ValueError:
                        file_posix = str(Path(exc["file"]).as_posix())
                    import_path_norm = re.sub(r"\s+", " ", exc["import_path"]).strip()
                    key = (source_id, exc["target_id"], file_posix, import_path_norm, False)
                    exception_keys.add(key)
        obsolete = set()
        for key in exception_keys:
            if key not in matched_exceptions:
                obsolete.add(key)
        assert len(obsolete) == 1, "Exception without observed violation should be detected as obsolete"
