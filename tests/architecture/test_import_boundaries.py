import ast
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

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

    def visit_If(self, node: ast.If):
        is_tc = False
        if isinstance(node.test, ast.Name) and node.test.id == "TYPE_CHECKING":
            is_tc = True
        elif isinstance(node.test, ast.Attribute):
            if (isinstance(node.test.value, ast.Name) and node.test.value.id == "typing" and
                    node.test.attr == "TYPE_CHECKING"):
                is_tc = True

        self.context_stack.append(self.in_type_checking)
        if is_tc:
            self.in_type_checking = True

        for stmt in node.body:
            self.visit(stmt)
        for stmt in node.orelse:
            self.visit(stmt)

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
        try:
            rel_path = fp.relative_to(src_chemuson)
            parts = list(rel_path.parts)
            if parts[-1] == "__init__.py":
                parts.pop()
            elif parts[-1].endswith(".py"):
                parts[-1] = parts[-1][:-3]
            if not parts or parts == ["__init__"]:
                return "chemuson"
            return "chemuson." + ".".join(parts)
        except ValueError:
            return ""

    def _resolve_import(self, imp_data: Dict[str, Any], fp: Path) -> Optional[tuple]:
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
            return (self._get_module_name(fp), best_match_id, imp_data["type_checking_only"])
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
            except SyntaxError:
                continue

            visitor = ImportVisitor(source, fp, self._get_module_name(fp))
            visitor.visit(tree)

            for imp_data in visitor.imports:
                resolved = self._resolve_import(imp_data, fp)
                if resolved:
                    source_id, target_id, tc_only = resolved
                    if target_id != m_id:
                        try:
                            node = imp_data["node"]
                            stmt = ast.get_source_segment(source, node)
                            if stmt is None:
                                stmt = ast.unparse(node)
                            stmt = re.sub(r"\s+", " ", stmt).strip()
                        except:
                            stmt = "unknown"

                        all_observed.append(ObservedImport(
                            source_id=source_id,
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
        formula_file = src_chemuson / "chemcalc" / "formula.py"

        imp_data = {
            "type": "from",
            "module": "valence",
            "level": 1,
            "names": ["implicit_h_count"],
            "lineno": 10,
            "node": ast.parse("from .valence import implicit_h_count").body[0],
            "type_checking_only": False
        }

        resolved = analyzer._resolve_import(imp_data, formula_file)
        assert resolved is not None
        source_id, target_id, tc_only = resolved
        assert source_id == "chemuson.chemcalc.formula"
        assert target_id == "M03"
        assert tc_only is False

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
        for imp in all_observed:
            assert imp.source_id != ""
            assert imp.target_id != ""
            assert imp.file.exists()
            assert imp.line > 0
            assert len(imp.statement) > 0
            assert imp.source_id != imp.target_id
            assert imp.source_id.startswith("chemuson")
            assert re.match(r"^M\d\d$", imp.target_id)
