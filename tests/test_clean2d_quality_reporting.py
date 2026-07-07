from __future__ import annotations

import pytest

from chemuson.clean2d import (
    Clean2DCandidate,
    Clean2DMode,
    Clean2DResult,
    rank_clean2d_candidates,
    summarize_clean2d_candidates,
)
from chemuson.clean2d.quality_reporting import (
    Clean2DQualityContractError,
    clean2d_quality_diagnostic,
    clean2d_result_quality_diagnostic,
    clean2d_safety_quality_diagnostic,
    normalize_clean2d_reporting_score,
    validate_clean2d_quality_diagnostic,
)
from chemuson.clean2d.safety import Clean2DQualityReport
from chemuson.core.model import MolGraph


def test_quality_diagnostic_accepts_applied_without_reason() -> None:
    diagnostic = clean2d_quality_diagnostic(
        state="applied",
        source="engine",
        score=0.85,
        internal={"raw_score": 12.0},
    )

    assert diagnostic == {
        "state": "applied",
        "score": 0.85,
        "source": "engine",
        "internal": {"raw_score": 12.0},
    }
    validate_clean2d_quality_diagnostic(diagnostic)


@pytest.mark.parametrize("state", ["rejected", "failed-controlled"])
def test_quality_diagnostic_requires_reason_for_rejected_and_failed_controlled(state: str) -> None:
    with pytest.raises(Clean2DQualityContractError):
        clean2d_quality_diagnostic(state=state, source="safety", score=0.0)


@pytest.mark.parametrize("state", ["no-op", "preserve-only"])
def test_quality_diagnostic_allows_optional_reason_for_noop_and_preserve_only(state: str) -> None:
    without_reason = clean2d_quality_diagnostic(state=state, source="engine", score=1.0)
    with_reason = clean2d_quality_diagnostic(
        state=state,
        source="engine",
        score=1.0,
        reason="worse-quality",
    )

    assert "reason" not in without_reason
    assert with_reason["reason"] == "worse-quality"


def test_quality_diagnostic_rejects_applied_reason_in_contract() -> None:
    with pytest.raises(Clean2DQualityContractError):
        clean2d_quality_diagnostic(
            state="applied",
            source="engine",
            score=1.0,
            reason="worse-quality",
        )


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (-10.0, 1.0),
        (0.0, 1.0),
        (50.0, 0.5),
        (100.0, 0.0),
        (250.0, 0.0),
    ],
)
def test_normalize_clean2d_reporting_score_is_bounded_and_lower_raw_is_better(
    raw: float,
    expected: float,
) -> None:
    assert normalize_clean2d_reporting_score(raw, lower_is_better=True, scale=100.0) == expected


def test_safety_report_adapter_maps_real_rejection_reason_and_preserves_internal_metadata() -> None:
    report = Clean2DQualityReport(
        atom_ids={1, 2},
        before={1: (0.0, 0.0), 2: (40.0, 0.0)},
        after={1: (0.0, 0.0), 2: (40.0, 0.0)},
        target_bond_length=40.0,
        passed=False,
        rejection_reason="colision_no_enlazada",
        min_nonbonded_after=2.0,
    )

    diagnostic = clean2d_safety_quality_diagnostic(report)

    assert diagnostic["state"] == "rejected"
    assert diagnostic["reason"] == "collision-risk"
    assert diagnostic["source"] == "safety"
    assert 0.0 <= diagnostic["score"] <= 1.0
    assert diagnostic["internal"]["rejection_reason"] == "colision_no_enlazada"
    assert diagnostic["internal"]["min_nonbonded_after"] == 2.0


def test_result_adapter_reports_candidate_source_score_state_and_internal_metadata() -> None:
    before = {1: (0.0, 0.0), 2: (40.0, 0.0)}
    candidate = Clean2DCandidate(
        source="rdkit_isolated",
        coords={1: (0.0, 1.0), 2: (40.0, 1.0)},
        score=12.5,
        novelty=1.0,
        metadata={"quality_class": "good", "visual_score": 12.5},
    )
    result = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords=before,
        candidates=(candidate,),
        selected=candidate,
    )

    diagnostic = clean2d_result_quality_diagnostic(result)

    assert diagnostic["state"] == "applied"
    assert diagnostic["source"] == "rdkit_isolated"
    assert "reason" not in diagnostic
    assert 0.0 <= diagnostic["score"] <= 1.0
    assert diagnostic["internal"]["candidate_sources"] == ("rdkit_isolated",)
    assert diagnostic["internal"]["selected_score"] == 12.5

    summary = summarize_clean2d_candidates(result)
    assert summary[0]["quality_diagnostic"]["state"] == "applied"
    assert summary[0]["quality_diagnostic"]["source"] == "rdkit_isolated"


def test_reporting_score_does_not_change_existing_candidate_selection() -> None:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    graph.add_atom("C", 40.0, 0.0, atom_id=2)
    graph.add_bond(1, 2, order=1)
    before = {1: (0.0, 0.0), 2: (40.0, 0.0)}
    candidates = [
        Clean2DCandidate(source="higher_raw", coords=before, score=100.0, novelty=0.0),
        Clean2DCandidate(source="lower_raw", coords=before, score=0.0, novelty=0.0),
    ]

    selected_before = rank_clean2d_candidates(graph, candidates, before, before.keys()).selected
    for candidate in candidates:
        clean2d_quality_diagnostic(
            state=candidate.outcome_state,
            source=candidate.source,
            score=normalize_clean2d_reporting_score(candidate.score, lower_is_better=True, scale=100.0),
            internal={"raw_score": candidate.score},
        )
    selected_after = rank_clean2d_candidates(graph, candidates, before, before.keys()).selected

    assert selected_before is not None
    assert selected_after is not None
    assert selected_before.source == selected_after.source
