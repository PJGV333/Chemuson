from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "docs" / "operational-resilience-audit.md"


def test_operational_resilience_ownership_is_explicit() -> None:
    text = AUDIT.read_text(encoding="utf-8").lower()
    for term in ("m22", "crash", "autosave", "m08", "m10", "m14", "telemetry"):
        assert term in text
    assert "no new module" in text
