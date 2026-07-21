import ast
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

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
                        try:
                            stmt = ast.get_source_segment(source, node)
                            if stmt is None:
                                stmt = ast.unparse(node)
                            stmt = re.sub(r"\s+", " ", stmt).strip()
                        except Exception:
                            stmt = ast.unparse(node)

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

    def test_type_checking_else_block_is_runtime(self):
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
        assert "chemuson.gui.canvas" in (
            tc_imports[0].get("module", "") or tc_imports[0].get("name", "")
        ), (
            f"Expected chemuson.gui.canvas in TYPE_CHECKING import, got {tc_imports[0]}"
        )

        assert len(rt_imports) == 1, (
            f"Expected exactly 1 runtime import in else block, found {len(rt_imports)}"
        )
        assert "MolGraph" in rt_imports[0]["names"][0], (
            f"Expected MolGraph in runtime import, got {rt_imports[0]['names']}"
        )
        assert "chemuson.core" in (
            rt_imports[0].get("module", "") or rt_imports[0].get("name", "")
        ), (
            f"Expected chemuson.core in runtime import, got {rt_imports[0]}"
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
