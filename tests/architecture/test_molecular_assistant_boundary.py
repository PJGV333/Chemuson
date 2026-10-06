from __future__ import annotations

import ast
import os
import subprocess
import sys
from dataclasses import fields
from pathlib import Path

import yaml

from chemuson.molecular_assistant import (
    MolecularAssistantRequest,
    MolecularTransformationRequest,
)


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "architecture" / "modules.yml"
SOURCE = ROOT / "src" / "chemuson" / "molecular_assistant"
CLEAN2D_SOURCE = ROOT / "src" / "chemuson" / "clean2d"


def _modules() -> dict[str, dict[str, object]]:
    data = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    return {module["id"]: module for module in data["modules"]}


def _import_names(tree: ast.AST) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def test_m23_catalog_has_only_core_and_chemio_as_dependencies():
    modules = _modules()
    assistant = modules["M23"]

    assert assistant["name"] == "molecular_assistant"
    assert assistant["paths"] == ["src/chemuson/molecular_assistant/"]
    assert assistant["current_dependencies"] == ["M00", "M01"]
    assert assistant["target_dependencies"] == ["M00", "M01"]
    assert set(assistant["forbidden_dependencies"]) >= {
        "M02",
        "M04",
        "M08",
        "M09",
        "M10",
        "M11",
        "M12",
        "M13",
        "M16",
        "M19",
    }
    assert assistant["temporary_exceptions"] == []
    assert assistant["circular_dependencies"] == []


def test_core_chemio_and_clean2d_forbid_reverse_ai_dependency():
    modules = _modules()
    for module_id in ("M00", "M01", "M02"):
        assert "M23" in modules[module_id]["forbidden_dependencies"]
        assert "M23" not in modules[module_id]["current_dependencies"]
        assert "M23" not in modules[module_id]["target_dependencies"]


def test_m23_source_imports_only_core_chemio_and_its_own_package():
    allowed = {"chemuson.core", "chemuson.chemio", "chemuson.molecular_assistant"}
    for path in SOURCE.rglob("*.py"):
        names = _import_names(ast.parse(path.read_text(encoding="utf-8")))
        chemuson_imports = {name for name in names if name == "chemuson" or name.startswith("chemuson.")}
        for imported in chemuson_imports:
            assert any(imported == prefix or imported.startswith(prefix + ".") for prefix in allowed), (
                f"Unexpected ChemUSON dependency in {path.relative_to(ROOT)}: {imported}"
            )
        assert not any(name == "tools" or name.startswith("tools.") for name in names)


def test_clean2d_source_has_no_ai_imports():
    for path in CLEAN2D_SOURCE.rglob("*.py"):
        names = _import_names(ast.parse(path.read_text(encoding="utf-8")))
        assert not any(
            name == "chemuson.molecular_assistant"
            or name.startswith("chemuson.molecular_assistant.")
            for name in names
        ), f"Clean2D imported the AI module in {path.relative_to(ROOT)}"


def test_ai_clean2d_orchestration_is_tool_only_and_does_not_add_reverse_imports():
    evaluator = ROOT / "tools" / "ai_clean2d_evaluation.py"
    imports = _import_names(ast.parse(evaluator.read_text(encoding="utf-8")))
    assert "chemuson.clean2d" in imports
    assert "chemuson.molecular_assistant" in imports
    assert not any(name.startswith("chemuson.gui") for name in imports)
    assert evaluator.parent == ROOT / "tools"
    test_clean2d_source_has_no_ai_imports()


def test_request_contains_no_canvas_document_or_mutation_context():
    assert [field.name for field in fields(MolecularAssistantRequest)] == ["description"]
    assert [field.name for field in fields(MolecularTransformationRequest)] == [
        "source_smiles",
        "instruction",
    ]


def test_package_import_does_not_load_gui_clean2d_chemname_or_rdkit():
    script = (
        "import sys; import chemuson.molecular_assistant; "
        "blocked = ('chemuson.gui', 'chemuson.clean2d', 'chemuson.chemname', 'rdkit'); "
        "loaded = [name for name in sys.modules if any(name == item or name.startswith(item + '.') "
        "for item in blocked)]; "
        "assert not loaded, loaded"
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run(
        [sys.executable, "-c", script],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
