from __future__ import annotations

import json

from tests.clean2d_regression.reports import build_baseline_report, canonicalize_baseline_report, compare_baseline_reports
from tests.clean2d_regression.review_policy import classify_baseline_diff, summarize_baseline_diff_review


def _single_field_diff(field: str, value, *, known_delicate: bool = False):
    left = canonicalize_baseline_report(build_baseline_report())
    right = canonicalize_baseline_report(left)
    index = 0
    if known_delicate:
        index = next(idx for idx, record in enumerate(right["records"]) if "known_delicate" in record["tags"])
    right["records"][index][field] = value
    return compare_baseline_reports(left, right), left


def test_empty_diff_produces_equivalent_review_summary() -> None:
    diff = compare_baseline_reports(build_baseline_report(), build_baseline_report())
    items = classify_baseline_diff(diff)
    summary = summarize_baseline_diff_review(items)

    assert items == ()
    assert summary["equivalent"] is True
    assert summary["item_count"] == 0
    json.dumps(summary, allow_nan=False, sort_keys=True)


def test_result_state_change_is_contract_risk() -> None:
    diff, report = _single_field_diff("result_state", "failed-controlled")
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "result_state"
    assert item.category == "contract"
    assert item.severity == "contract-risk"


def test_stable_reason_change_needs_review() -> None:
    diff, report = _single_field_diff("stable_reason", "backend-failure")
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "stable_reason"
    assert item.category == "contract"
    assert item.severity == "needs-review"


def test_selected_source_change_is_routing_review() -> None:
    diff, report = _single_field_diff("selected_source", "artificial-source")
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "selected_source"
    assert item.category == "routing"
    assert item.severity == "needs-review"


def test_candidate_sources_change_is_routing_review() -> None:
    diff, report = _single_field_diff("candidate_sources", ["artificial-source"])
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "candidate_sources"
    assert item.category == "routing"
    assert item.severity == "needs-review"


def test_metric_change_is_geometry_diagnostic_observational_only() -> None:
    left = canonicalize_baseline_report(build_baseline_report())
    right = canonicalize_baseline_report(left)
    right["records"][0]["metrics"]["before"]["length_rms_error"] += 1.0
    diff = compare_baseline_reports(left, right)
    item = classify_baseline_diff(diff, report=left)[0]

    assert item.field == "metrics"
    assert item.category == "geometry-diagnostic"
    assert item.severity == "geometry-risk"
    assert item.observational_only is True


def test_snapshot_change_is_diagnostic_review() -> None:
    diff, report = _single_field_diff("snapshot", {"changed": True})
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "snapshot"
    assert item.category == "diagnostic"
    assert item.severity == "needs-review"


def test_policy_evidence_change_is_complex_policy_review() -> None:
    diff, report = _single_field_diff("policy_evidence", {"changed": True})
    item = classify_baseline_diff(diff, report=report)[0]

    assert item.field == "policy_evidence"
    assert item.category == "complex-policy"
    assert item.severity == "needs-review"


def test_known_delicate_change_remains_visible_and_labeled() -> None:
    diff, report = _single_field_diff("stable_reason", "backend-failure", known_delicate=True)
    item = classify_baseline_diff(diff, report=report)[0]
    summary = summarize_baseline_diff_review((item,))

    assert item.known_delicate_related is True
    assert item.severity == "known-delicate-change"
    assert summary["known_delicate_related_count"] == 1
    assert summary["item_count"] == 1


def test_review_summary_is_json_stable() -> None:
    diff, report = _single_field_diff("selected_source", "artificial-source")
    items = classify_baseline_diff(diff, report=report)
    summary = summarize_baseline_diff_review(items)

    json.dumps(summary, allow_nan=False, sort_keys=True)
    assert summary["observational_only"] is True
    assert summary["equivalent"] is False


def test_review_policy_does_not_modify_reports() -> None:
    report = canonicalize_baseline_report(build_baseline_report())
    before = json.dumps(report, allow_nan=False, sort_keys=True)
    diff, _ = _single_field_diff("selected_source", "artificial-source")

    classify_baseline_diff(diff, report=report)

    assert json.dumps(report, allow_nan=False, sort_keys=True) == before
