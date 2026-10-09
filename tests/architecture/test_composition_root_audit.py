from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "docs" / "composition-root-audit.md"


def test_composition_root_audit_records_consolidated_m19() -> None:
    text = AUDIT.read_text(encoding="utf-8")
    assert "M19" in text
    assert "consolidated" in text.lower()
    assert "no structural change required" in text
    assert "M24" in text


def test_bootstrap_catalog_tracks_private_packaging_smoke_dependencies() -> None:
    catalog = yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text())
    bootstrap = next(module for module in catalog["modules"] if module["id"] == "M19")
    assert bootstrap["paths"] == ["src/chemuson/__main__.py", "src/chemuson/app/"]
    expected_dependencies = {"M01", "M04", "M08", "M18", "M22"}
    assert set(bootstrap["current_dependencies"]) == expected_dependencies
    assert set(bootstrap["target_dependencies"]) == expected_dependencies
    assert "privado" in bootstrap["notes"]
    assert "M04" in bootstrap["notes"]
