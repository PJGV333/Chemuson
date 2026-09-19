from __future__ import annotations

import ast
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "src" / "chemuson" / "resilience" / "recovery.py"
CONTROLLER = ROOT / "src" / "chemuson" / "gui" / "controllers" / "recovery_controller.py"


def _module(module_id: str) -> dict:
    catalog = yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text(encoding="utf-8"))
    return next(module for module in catalog["modules"] if module["id"] == module_id)


def test_recovery_policy_has_one_canonical_owner() -> None:
    assert RECOVERY.exists()
    tree = ast.parse(RECOVERY.read_text(encoding="utf-8"))
    names = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert names == {
        "read_autosave_metadata",
        "list_autosave_entries",
        "archive_autosave",
    }

    controller_text = CONTROLLER.read_text(encoding="utf-8")
    assert "chemuson.resilience.recovery" in controller_text
    assert "open(" not in controller_text
    assert "os.listdir" not in controller_text
    assert "os.replace" not in controller_text
    assert "json.load" not in controller_text
    assert "datetime.now" not in controller_text


def test_recovery_catalog_and_dependencies_are_explicit() -> None:
    resilience = _module("M22")
    controllers = _module("M10")
    utils = _module("M15")
    assert "recovery" in resilience["internal_api"]
    assert "read_autosave_metadata" in resilience["public_api"]
    assert "list_autosave_entries" in resilience["public_api"]
    assert "archive_autosave" in resilience["public_api"]
    assert "tests/test_recovery_policy.py" in resilience["tests"]
    assert "tests/architecture/test_recovery_boundary.py" in resilience["tests"]
    assert "M22" in controllers["current_dependencies"]
    assert "M22" in controllers["target_dependencies"]
    assert resilience["current_dependencies"] == []
    assert resilience["target_dependencies"] == []
    assert resilience["temporary_exceptions"] == []
    assert resilience["circular_dependencies"] == []
    assert "Shims históricos" in utils["responsibility"]
    assert "autosave" in utils["notes"]
