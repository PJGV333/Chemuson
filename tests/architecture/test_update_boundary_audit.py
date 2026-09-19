from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "docs" / "update-boundary-audit.md"


def test_update_boundary_audit_records_no_extraction() -> None:
    text = AUDIT.read_text(encoding="utf-8")
    assert "M14" in text
    assert "audited / no structural change required" in text
    assert "M23" in text
    assert "no new module" in text.lower()


def test_update_catalog_remains_m14() -> None:
    catalog = yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text())
    update = next(module for module in catalog["modules"] if module["id"] == "M14")
    assert update["name"] == "update"
    assert update["paths"] == ["src/chemuson/update/"]
