"""Tests arquitectónicos: validación del catálogo de módulos (Phase 1).

Valida la estructura de `architecture/modules.yml` sin importar módulos de producción.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"

MODULE_ID_PATTERN = re.compile(r"^M\d\d$")
EXPECTED_IDS = {f"M{i:02d}" for i in range(20)}
VALID_STATUSES = {"stable", "evolving", "legacy"}
VALID_RISK_LEVELS = {"low", "medium", "high"}
VALID_SEVERITIES = {"low", "medium", "high"}
REQUIRED_FIELDS = [
    "id",
    "name",
    "title",
    "responsibility",
    "paths",
    "status",
    "risk_level",
    "current_dependencies",
    "target_dependencies",
    "forbidden_dependencies",
    "temporary_exceptions",
    "circular_dependencies",
    "entry_points",
    "tests",
    "public_api",
    "internal_api",
]
REQUIRED_EXCEPTION_FIELDS = {
    "source_id",
    "target_id",
    "file",
    "import_path",
    "reason",
    "debt_ref",
    "elimination_condition",
    "type_checking_only",
}


# ---------------------------------------------------------------------------
# Fixtures / helpers
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


# ---------------------------------------------------------------------------
# 1. Carga del catálogo
# ---------------------------------------------------------------------------

class TestCatalogLoad:
    """Validate that the catalog file exists and is parseable."""

    def test_catalog_file_exists(self):
        assert CATALOG_PATH.exists(), f"Catalog file not found: {CATALOG_PATH}"

    def test_yaml_is_valid(self, catalog):
        assert isinstance(catalog, dict), "YAML root is not a mapping"

    def test_root_contains_modules_key(self, catalog):
        assert "modules" in catalog, "Root mapping must contain 'modules'"

    def test_modules_is_list(self, catalog):
        assert isinstance(catalog["modules"], list), "'modules' must be a list"


# ---------------------------------------------------------------------------
# 2. IDs
# ---------------------------------------------------------------------------

class TestModuleIds:
    """Validate module identification: exactly 20 modules, M00-M19, unique, correct format."""

    def test_exactly_20_modules(self, modules):
        assert len(modules) == 20, f"Expected 20 modules, got {len(modules)}"

    def test_ids_are_m00_to_m19(self, modules):
        ids = {m["id"] for m in modules}
        assert ids == EXPECTED_IDS, f"IDs mismatch: {ids}"

    def test_no_duplicate_ids(self, modules):
        ids = [m["id"] for m in modules]
        assert len(ids) == len(set(ids)), "Duplicate module IDs found"

    def test_each_id_matches_pattern(self, modules):
        for m in modules:
            assert MODULE_ID_PATTERN.match(m["id"]), (
                f"Module ID '{m['id']}' does not match pattern ^M\\d\\d$"
            )


# ---------------------------------------------------------------------------
# 3. Campos obligatorios
# ---------------------------------------------------------------------------

class TestMandatoryFields:
    """Validate that every module has all required fields with correct types."""

    def test_all_required_fields_present(self, modules):
        for m in modules:
            for field in REQUIRED_FIELDS:
                assert field in m, (
                    f"Module '{m['id']}' missing required field '{field}'"
                )

    def test_string_fields_are_strings(self, modules):
        string_fields = {"id", "name", "title", "responsibility", "status", "risk_level"}
        for m in modules:
            for field in string_fields:
                assert isinstance(m[field], str), (
                    f"Module '{m['id']}', field '{field}' must be str, got {type(m[field]).__name__}"
                )

    def test_list_fields_are_lists(self, modules):
        list_fields = [
            "paths",
            "current_dependencies",
            "target_dependencies",
            "forbidden_dependencies",
            "temporary_exceptions",
            "circular_dependencies",
            "entry_points",
            "tests",
            "public_api",
            "internal_api",
        ]
        for m in modules:
            for field in list_fields:
                assert isinstance(m[field], list), (
                    f"Module '{m['id']}', field '{field}' must be list, got {type(m[field]).__name__}"
                )

    def test_status_values_valid(self, modules):
        for m in modules:
            assert m["status"] in VALID_STATUSES, (
                f"Module '{m['id']}', status '{m['status']}' not in {VALID_STATUSES}"
            )

    def test_risk_level_values_valid(self, modules):
        for m in modules:
            assert m["risk_level"] in VALID_RISK_LEVELS, (
                f"Module '{m['id']}', risk_level '{m['risk_level']}' not in {VALID_RISK_LEVELS}"
            )


# ---------------------------------------------------------------------------
# 4. Rutas
# ---------------------------------------------------------------------------

class TestModulePaths:
    """Validate that declared paths exist and no path overlaps."""

    def _resolve_path(self, path_str):
        """Resolve a path relative to the repository root."""
        return (REPO_ROOT / path_str).resolve()

    def test_no_absolute_paths(self, modules):
        for m in modules:
            for p in m["paths"]:
                assert not Path(p).is_absolute(), (
                    f"Module '{m['id']}', path '{p}' must not be absolute"
                )

    def test_no_path_contains_dotdot(self, modules):
        for m in modules:
            for p in m["paths"]:
                assert ".." not in p, (
                    f"Module '{m['id']}', path '{p}' must not contain '..'"
                )

    def test_all_paths_exist(self, modules):
        for m in modules:
            for p in m["paths"]:
                resolved = self._resolve_path(p)
                assert resolved.exists(), (
                    f"Module '{m['id']}', path '{p}' does not exist"
                )

    def test_python_files_belong_to_at_most_one_module(self, modules):
        """Each .py file in src/chemuson/ belongs to at most one module."""
        module_paths = {}  # module_id -> set of resolved paths
        for m in modules:
            module_paths[m["id"]] = set()
            for p in m["paths"]:
                resolved = self._resolve_path(p)
                if resolved.is_dir():
                    for py_file in resolved.rglob("*.py"):
                        if "__pycache__" in py_file.parts or py_file.suffix == ".pyc":
                            continue
                        module_paths[m["id"]].add(py_file.resolve())
                else:
                    module_paths[m["id"]].add(resolved)

        # Check for overlaps
        all_files = {}  # resolved_path -> list of module_ids
        for mid, paths in module_paths.items():
            for fp in paths:
                if fp not in all_files:
                    all_files[fp] = []
                all_files[fp].append(mid)

        for fp, ids in all_files.items():
            if len(ids) > 1:
                assert False, (
                    f"File '{fp}' belongs to multiple modules: {ids}"
                )


# ---------------------------------------------------------------------------
# 5. Dependencias
# ---------------------------------------------------------------------------

class TestModuleDependencies:
    """Validate dependency fields: valid IDs, no duplicates, no self-deps, no overlaps."""

    def test_current_dependencies_valid_ids(self, modules):
        for m in modules:
            for dep_id in m["current_dependencies"]:
                assert dep_id in EXPECTED_IDS, (
                    f"Module '{m['id']}', 'current_dependencies' contains "
                    f"uncatalogued ID '{dep_id}'"
                )

    def test_target_dependencies_valid_ids(self, modules):
        for m in modules:
            for dep_id in m["target_dependencies"]:
                assert dep_id in EXPECTED_IDS, (
                    f"Module '{m['id']}', 'target_dependencies' contains "
                    f"uncatalogued ID '{dep_id}'"
                )

    def test_forbidden_dependencies_valid_ids(self, modules):
        for m in modules:
            for dep_id in m["forbidden_dependencies"]:
                assert dep_id in EXPECTED_IDS, (
                    f"Module '{m['id']}', 'forbidden_dependencies' contains "
                    f"uncatalogued ID '{dep_id}'"
                )

    def test_no_duplicate_ids_in_current_dependencies(self, modules):
        for m in modules:
            ids = m["current_dependencies"]
            assert len(ids) == len(set(ids)), (
                f"Module '{m['id']}', current_dependencies has duplicates"
            )

    def test_no_duplicate_ids_in_target_dependencies(self, modules):
        for m in modules:
            ids = m["target_dependencies"]
            assert len(ids) == len(set(ids)), (
                f"Module '{m['id']}', target_dependencies has duplicates"
            )

    def test_no_duplicate_ids_in_forbidden_dependencies(self, modules):
        for m in modules:
            ids = m["forbidden_dependencies"]
            assert len(ids) == len(set(ids)), (
                f"Module '{m['id']}', forbidden_dependencies has duplicates"
            )

    def test_no_self_dependency(self, modules):
        """Identify explicitly which dependency list contains self-references."""
        for m in modules:
            if m["id"] in m["current_dependencies"]:
                assert False, (
                    f"Module '{m['id']}' has self-dependency in "
                    f"current_dependencies"
                )
            if m["id"] in m["target_dependencies"]:
                assert False, (
                    f"Module '{m['id']}' has self-dependency in "
                    f"target_dependencies"
                )

    def test_target_and_forbidden_no_overlap(self, modules):
        for m in modules:
            target = set(m["target_dependencies"])
            forbidden = set(m["forbidden_dependencies"])
            overlap = target & forbidden
            assert not overlap, (
                f"Module '{m['id']}', target and forbidden dependencies overlap: {overlap}"
            )

    def test_m03_depends_only_on_m00_without_exceptions(self, modules):
        chemcalc = next(module for module in modules if module["id"] == "M03")
        assert chemcalc["current_dependencies"] == ["M00"]
        assert chemcalc["temporary_exceptions"] == []

    def test_m01_m02_dependency_direction_and_ownership(self, modules):
        m01 = next(module for module in modules if module["id"] == "M01")
        m02 = next(module for module in modules if module["id"] == "M02")

        assert m01["current_dependencies"] == ["M00"]
        assert m01["target_dependencies"] == ["M00"]
        assert "M02" not in m01["current_dependencies"]
        assert "M02" not in m01["target_dependencies"]
        assert m01["circular_dependencies"] == []
        assert len(m01["temporary_exceptions"]) == 1
        exception = m01["temporary_exceptions"][0]
        assert exception["source_id"] == "M01"
        assert exception["target_id"] == "M09"
        assert exception["type_checking_only"] is True
        assert exception["file"] == "src/chemuson/chemio/persistence.py"
        assert "depiction_candidates" not in m01["internal_api"]
        assert "depiction_quality" not in m01["internal_api"]
        assert "imported_depiction" not in m01["internal_api"]

        assert set(m02["current_dependencies"]) == {"M00", "M01"}
        assert set(m02["target_dependencies"]) == {"M00", "M01"}
        assert m02["circular_dependencies"] == []
        assert m02["temporary_exceptions"] == []
        assert "depiction_quality" in m02["internal_api"]
        assert "imported_depiction" in m02["internal_api"]
        for symbol in (
            "DepictionCandidate",
            "smiles_to_depiction_candidates",
            "smiles_to_molgraph_best_depiction",
            "smiles_to_molgraph_best_depiction_with_report",
        ):
            assert symbol in m02["public_api"]
        for test_path in (
            "tests/test_smiles_depiction_candidates.py",
            "tests/test_block_unwrap_depiction.py",
            "tests/test_scaffold_depiction.py",
            "tests/test_smiles_stereo_import.py",
        ):
            assert test_path in m02["tests"]

    def test_m00_registers_molecular_view_coverage(self, modules):
        core = next(module for module in modules if module["id"] == "M00")
        assert "tests/test_molecular_view.py" in core["tests"]

    def test_m00_to_m18_no_m19_dependency(self, modules):
        """Modules M00-M18 must not reference M19 in current or target dependencies."""
        for m in modules:
            if m["id"] == "M19":
                continue
            for dep_list, list_name in [
                (m["current_dependencies"], "current_dependencies"),
                (m["target_dependencies"], "target_dependencies"),
            ]:
                if "M19" in dep_list:
                    assert False, (
                        f"Module '{m['id']}' references M19 in "
                        f"{list_name}"
                    )


# ---------------------------------------------------------------------------
# 6. Excepciones
# ---------------------------------------------------------------------------

class TestTemporaryExceptions:
    """Validate the structure of temporary_exceptions entries."""

    def test_exception_has_required_fields(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                missing = REQUIRED_EXCEPTION_FIELDS - set(exc.keys())
                assert not missing, (
                    f"Module '{m['id']}', exception missing fields: {missing}"
                )

    def test_exception_source_id_is_valid(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                assert exc["source_id"] in EXPECTED_IDS, (
                    f"Module '{m['id']}', exception source_id '{exc['source_id']}' "
                    f"is not in catalog"
                )

    def test_exception_target_id_is_valid(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                assert exc["target_id"] in EXPECTED_IDS, (
                    f"Module '{m['id']}', exception target_id '{exc['target_id']}' "
                    f"is not in catalog"
                )

    def test_exception_source_id_matches_module(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                assert exc["source_id"] == m["id"], (
                    f"Module '{m['id']}', exception source_id '{exc['source_id']}' "
                    f"does not match module ID"
                )

    def test_exception_file_exists(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                resolved = (REPO_ROOT / exc["file"]).resolve()
                assert resolved.exists(), (
                    f"Module '{m['id']}', exception file '{exc['file']}' does not exist"
                )

    def test_exception_string_fields_non_empty(self, modules):
        string_fields = {"import_path", "reason", "debt_ref", "elimination_condition"}
        for m in modules:
            for exc in m["temporary_exceptions"]:
                for field in string_fields:
                    assert isinstance(exc[field], str), (
                        f"Module '{m['id']}', exception '{field}' must be str"
                    )
                    assert exc[field], (
                        f"Module '{m['id']}', exception '{field}' must not be empty"
                    )

    def test_exception_type_checking_only_is_bool(self, modules):
        for m in modules:
            for exc in m["temporary_exceptions"]:
                assert isinstance(exc["type_checking_only"], bool), (
                    f"Module '{m['id']}', exception 'type_checking_only' must be bool"
                )

    def test_no_duplicate_exceptions(self, modules):
        """No two exceptions share the same (source_id, target_id, file, import_path, type_checking_only)."""
        seen = set()
        for m in modules:
            for exc in m["temporary_exceptions"]:
                key = (
                    exc["source_id"],
                    exc["target_id"],
                    exc["file"],
                    exc["import_path"],
                    exc["type_checking_only"],
                )
                assert key not in seen, (
                    f"Duplicate exception: {key}"
                )
                seen.add(key)


# ---------------------------------------------------------------------------
# 7. Ciclos
# ---------------------------------------------------------------------------

class TestCircularDependencies:
    """Validate the structure of circular_dependencies entries."""

    def test_no_circular_dependencies_declared(self, modules):
        """The final catalog declares no circular dependencies."""
        all_cycles = [cycle for module in modules for cycle in module["circular_dependencies"]]
        assert all_cycles == []

    def test_cycle_severity_valid(self, modules):
        for m in modules:
            for cycle in m["circular_dependencies"]:
                assert cycle["severity"] in VALID_SEVERITIES, (
                    f"Cycle severity '{cycle['severity']}' not in {VALID_SEVERITIES}"
                )

    def test_cycle_has_required_fields(self, modules):
        for m in modules:
            for cycle in m["circular_dependencies"]:
                for field in ("modules", "edges", "severity", "resolution_plan"):
                    assert field in cycle, (
                        f"Cycle in module '{m['id']}' missing field '{field}'"
                    )

    def test_cycle_edges_have_required_fields(self, modules):
        for m in modules:
            for cycle in m["circular_dependencies"]:
                for edge in cycle["edges"]:
                    for field in ("source", "target", "file", "import_path"):
                        assert field in edge, (
                            f"Cycle edge in module '{m['id']}' missing field '{field}'"
                        )

    def test_cycle_edge_ids_valid(self, modules):
        for m in modules:
            for cycle in m["circular_dependencies"]:
                for edge in cycle["edges"]:
                    assert edge["source"] in EXPECTED_IDS, (
                        f"Cycle edge in module '{m['id']}', "
                        f"invalid source ID '{edge['source']}'"
                    )
                    assert edge["target"] in EXPECTED_IDS, (
                        f"Cycle edge in module '{m['id']}', "
                        f"invalid target ID '{edge['target']}'"
                    )

    def test_cycle_modules_are_valid_ids(self, modules):
        """Each module ID in a cycle must belong to EXPECTED_IDS."""
        for m in modules:
            for cycle in m["circular_dependencies"]:
                for mod_id in cycle["modules"]:
                    assert mod_id in EXPECTED_IDS, (
                        f"Cycle in module '{m['id']}', cycle module '{mod_id}' "
                        f"is not in catalog"
                    )

    def test_cycle_edge_ids_in_cycle_modules(self, modules):
        """Each edge source/target must be in the cycle's module set."""
        for m in modules:
            for cycle in m["circular_dependencies"]:
                module_set = cycle["modules"]
                for edge in cycle["edges"]:
                    assert edge["source"] in module_set, (
                        f"Cycle in module '{m['id']}', edge source '{edge['source']}' "
                        f"not in cycle modules {module_set}"
                    )
                    assert edge["target"] in module_set, (
                        f"Cycle in module '{m['id']}', edge target '{edge['target']}' "
                        f"not in cycle modules {module_set}"
                    )

    def test_cycle_edge_paths_exist(self, modules):
        for m in modules:
            for cycle in m["circular_dependencies"]:
                for edge in cycle["edges"]:
                    resolved = (REPO_ROOT / edge["file"]).resolve()
                    assert resolved.exists(), (
                        f"Cycle edge in module '{m['id']}', path '{edge['file']}' "
                        f"does not exist"
                    )

    def test_cycle_appears_once(self, modules):
        """Each cycle should appear exactly once."""
        all_cycles = []
        for m in modules:
            all_cycles.extend(m["circular_dependencies"])

        cycle_keys = []
        for cycle in all_cycles:
            edges_sorted = tuple(
                tuple(sorted(e.items())) for e in sorted(cycle["edges"], key=lambda e: e["file"])
            )
            key = (
                frozenset(cycle["modules"]),
                cycle["severity"],
                edges_sorted,
            )
            cycle_keys.append(key)

        assert len(cycle_keys) == len(set(cycle_keys)), (
            "Duplicate cycle found"
        )
