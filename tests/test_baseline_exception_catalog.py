"""Fail-closed validation for the release baseline exception catalog."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "packaging/release/validate_baseline_exceptions.py"
SPEC = importlib.util.spec_from_file_location("validate_baseline_exceptions", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_versioned_baseline_exception_catalog_is_exact_and_auditable() -> None:
    path = ROOT / "docs/release/KNOWN_BASELINE_EXCEPTIONS.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert MODULE.validate_catalog(data, ROOT) == 4


def test_catalog_rejects_wildcards_and_unknown_cases() -> None:
    path = ROOT / "docs/release/KNOWN_BASELINE_EXCEPTIONS.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["exceptions"][0]["case"] = "tests/test_*.py::test_anything"

    with pytest.raises(ValueError, match="wildcards"):
        MODULE.validate_catalog(data, ROOT)


def test_catalog_rejects_missing_evidence_and_blanket_skip_policy() -> None:
    path = ROOT / "docs/release/KNOWN_BASELINE_EXCEPTIONS.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["exceptions"][0]["evidence"] = "missing/evidence.md"

    with pytest.raises(ValueError, match="existing repository file"):
        MODULE.validate_catalog(data, ROOT)

    data = json.loads(path.read_text(encoding="utf-8"))
    data["policy"]["automatic_skips"] = True
    with pytest.raises(ValueError, match="fail-closed"):
        MODULE.validate_catalog(data, ROOT)
