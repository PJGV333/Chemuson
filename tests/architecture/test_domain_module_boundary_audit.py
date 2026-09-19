from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
AUDIT_PATH = ROOT / "docs" / "module-boundary-audit.md"
AUDITED_IDS = ("M05", "M06", "M07", "M14", "M16", "M17", "M18", "M19", "M20")


def _catalog() -> dict:
    return yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text())


def test_existing_domain_audit_records_every_audited_module() -> None:
    text = AUDIT_PATH.read_text()
    for module_id in AUDITED_IDS:
        assert f"| {module_id} |" in text
        assert f"{module_id}: audited / no structural change required" in text


def test_existing_domain_audit_matches_zero_debt_catalog_state() -> None:
    modules = {entry["id"]: entry for entry in _catalog()["modules"]}
    for module_id in AUDITED_IDS:
        entry = modules[module_id]
        assert entry["temporary_exceptions"] == []
        assert entry["circular_dependencies"] == []
        assert entry["current_dependencies"] == entry["target_dependencies"]
