"""Acceptance tests del runner de nomenclatura ChemName."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

from tools.chemname_acceptance import (  # noqa: E402
    DEFAULT_CASES_PATH,
    evaluate_case,
    is_rdkit_worker_available,
    load_cases,
    run_acceptance,
)

CASES = load_cases(DEFAULT_CASES_PATH)
CASE_IDS = [str(case.get("id", "")) for case in CASES]
RDKIT_WORKER_OK = is_rdkit_worker_available()


def test_acceptance_dataset_size_and_shape() -> None:
    """El dataset de acceptance debe mantenerse entre 50 y 100 casos."""
    assert 50 <= len(CASES) <= 100
    assert all(case_id.strip() for case_id in CASE_IDS)


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_acceptance_case(case: dict) -> None:
    """Cada caso debe pasar o marcarse como skipped cuando aplique."""
    if bool(case.get("skip_if_no_rdkit", False)) and not RDKIT_WORKER_OK:
        pytest.skip("RDKit worker no disponible")

    result = evaluate_case(case, DEFAULT_CASES_PATH, rdkit_available=RDKIT_WORKER_OK)
    if result.get("status") == "skip":
        pytest.skip(str(result.get("reason", "skipped")))
    assert result.get("status") == "pass", json.dumps(result, ensure_ascii=False, indent=2)


def test_acceptance_runner_writes_reports(tmp_path: Path) -> None:
    """El runner debe escribir reporte JSON y JUnit XML."""
    report_path = tmp_path / "chemname_acceptance_report.json"
    junit_path = tmp_path / "chemname_acceptance_report.xml"

    summary = run_acceptance(
        cases_path=DEFAULT_CASES_PATH,
        report_path=report_path,
        junit_path=junit_path,
        top_slowest=5,
        quiet=True,
    )

    assert report_path.exists()
    assert junit_path.exists()
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert int(payload.get("total", 0)) == len(CASES)
    assert int(summary.get("failed", 0)) == 0
    assert int(summary.get("errors", 0)) == 0

