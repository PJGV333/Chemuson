from __future__ import annotations

import copy
import inspect
import json

from tests.clean2d_regression import blocking_policy
from tests.clean2d_regression.blocking_policy import (
    ALLOWED_KNOWN_DELICATE,
    BLOCKING_CANDIDATE,
    CLEAN,
    REVIEW_REQUIRED,
    evaluate_diff_blocking_policy,
)


def _summary(*items: dict) -> dict:
    return {
        "equivalent": not items,
        "observational_only": True,
        "item_count": len(items),
        "items": list(items),
    }


def _item(severity: str, **overrides) -> dict:
    item = {
        "case_name": "ethanol",
        "field": "result_state",
        "category": "contract",
        "severity": severity,
        "known_delicate_related": False,
        "observational_only": True,
        "message": f"changed: {severity}",
    }
    item.update(overrides)
    return item


def test_empty_review_summary_produces_clean() -> None:
    decision = evaluate_diff_blocking_policy(_summary())

    assert decision.decision == CLEAN
    assert decision.blocking_candidate_count == 0
    assert decision.review_required_count == 0
    assert decision.informational_count == 0
    assert decision.reasons == ()
    assert decision.observational_only is True


def test_contract_risk_produces_blocking_candidate() -> None:
    decision = evaluate_diff_blocking_policy(_summary(_item("contract-risk")))

    assert decision.decision == BLOCKING_CANDIDATE
    assert decision.blocking_candidate_count == 1
    assert decision.reasons[0]["policy_classification"] == BLOCKING_CANDIDATE


def test_needs_review_produces_review_required() -> None:
    decision = evaluate_diff_blocking_policy(
        _summary(_item("needs-review", field="selected_source", category="routing"))
    )

    assert decision.decision == REVIEW_REQUIRED
    assert decision.review_required_count == 1


def test_geometry_risk_produces_review_required_and_observational_only() -> None:
    decision = evaluate_diff_blocking_policy(
        _summary(_item("geometry-risk", field="metrics", category="geometry-diagnostic"))
    )

    assert decision.decision == REVIEW_REQUIRED
    assert decision.review_required_count == 1
    assert decision.observational_only is True
    assert decision.reasons[0]["observational_only"] is True


def test_only_known_delicate_change_produces_allowed_known_delicate() -> None:
    decision = evaluate_diff_blocking_policy(
        _summary(_item("known-delicate-change", known_delicate_related=True))
    )

    assert decision.decision == ALLOWED_KNOWN_DELICATE
    assert decision.allowed_known_delicate_count == 1
    assert decision.reasons[0]["known_delicate_related"] is True


def test_known_delicate_change_does_not_hide_contract_risk() -> None:
    decision = evaluate_diff_blocking_policy(
        _summary(
            _item("known-delicate-change", case_name="delicate", known_delicate_related=True),
            _item("contract-risk", case_name="benzene"),
        )
    )

    assert decision.decision == BLOCKING_CANDIDATE
    assert decision.allowed_known_delicate_count == 1
    assert decision.blocking_candidate_count == 1
    assert [reason["policy_classification"] for reason in decision.reasons] == [
        ALLOWED_KNOWN_DELICATE,
        BLOCKING_CANDIDATE,
    ]


def test_decision_serializes_with_json_stability_options() -> None:
    decision = evaluate_diff_blocking_policy(_summary(_item("contract-risk")))

    encoded = json.dumps(decision.as_dict(), allow_nan=False, sort_keys=True)

    assert json.loads(encoded)["decision"] == BLOCKING_CANDIDATE


def test_policy_does_not_mutate_input_review_summary() -> None:
    review_summary = _summary(_item("needs-review", field="snapshot", category="diagnostic"))
    before = copy.deepcopy(review_summary)

    evaluate_diff_blocking_policy(review_summary)

    assert review_summary == before


def test_policy_is_test_only_and_does_not_import_production_modules() -> None:
    source = inspect.getsource(blocking_policy)

    assert "/tests/clean2d_regression/blocking_policy.py" in blocking_policy.__file__.replace("\\", "/")
    assert "from chemuson" not in source
    assert "import chemuson" not in source
    assert "from app" not in source
    assert "import app" not in source


def test_policy_does_not_introduce_automatic_gate() -> None:
    source = inspect.getsource(blocking_policy)
    decision = evaluate_diff_blocking_policy(_summary(_item("contract-risk")))

    assert decision.decision == BLOCKING_CANDIDATE
    assert decision.observational_only is True
    assert "pytest.fail" not in source
    assert "SystemExit" not in source
    assert "subprocess" not in source
