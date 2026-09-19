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


def test_bootstrap_catalog_ownership_is_unchanged() -> None:
    catalog = yaml.safe_load((ROOT / "architecture" / "modules.yml").read_text())
    bootstrap = next(module for module in catalog["modules"] if module["id"] == "M19")
    assert bootstrap["paths"] == ["src/chemuson/__main__.py", "src/chemuson/app/"]
