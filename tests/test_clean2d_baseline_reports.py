from __future__ import annotations

import json

from tests.clean2d_regression.cases import get_regression_cases
from tests.clean2d_regression.metrics import metric_tolerance
from tests.clean2d_regression.reports import (
    CLEAN2D_BASELINE_REPORT_SCHEMA,
    CLEAN2D_BASELINE_REPORT_VERSION,
    baseline_report_to_json,
    build_baseline_report,
    canonicalize_baseline_report,
    compare_baseline_reports,
)


def test_baseline_report_is_json_serializable() -> None:
    report = build_baseline_report(metadata={"run": "unit"})
    data = canonicalize_baseline_report(report)

    assert data["schema"] == CLEAN2D_BASELINE_REPORT_SCHEMA
    assert data["version"] == CLEAN2D_BASELINE_REPORT_VERSION
    assert data["metadata"] == {"run": "unit"}
    assert json.loads(baseline_report_to_json(report)) == data


def test_report_contains_one_record_per_corpus_case_sorted_by_name() -> None:
    report = canonicalize_baseline_report(build_baseline_report())
    names = [record["case_name"] for record in report["records"]]

    assert len(names) == len(get_regression_cases())
    assert names == sorted(names)
    assert set(names) == {case.name for case in get_regression_cases()}


def test_report_summary_counts_states_families_and_tags() -> None:
    report = canonicalize_baseline_report(build_baseline_report())
    records = report["records"]
    summary = report["summary"]

    assert summary["case_count"] == len(records)
    assert sum(summary["result_states"].values()) == len(records)
    assert sum(summary["families"].values()) == len(records)
    assert summary["tags"]["known_delicate"] == sum("known_delicate" in record["tags"] for record in records)
    assert summary["tags"]["complex_policy_guard"] == sum("complex_policy_guard" in record["tags"] for record in records)
    assert summary["tags"]["selection_boundary"] == sum("selection_boundary" in record["tags"] for record in records)
    assert summary["tags"]["stereo_sensitive"] == sum("stereo_sensitive" in record["tags"] for record in records)
    assert summary["tags"]["charged"] == sum("charged" in record["tags"] for record in records)


def test_consecutive_reports_are_equivalent() -> None:
    left = build_baseline_report()
    right = build_baseline_report()

    diff = compare_baseline_reports(left, right)

    assert diff == {
        "equivalent": True,
        "observational_only": True,
        "added_cases": [],
        "removed_cases": [],
        "changed_cases": [],
    }


def test_report_diff_detects_observable_field_changes() -> None:
    left = canonicalize_baseline_report(build_baseline_report())
    right = canonicalize_baseline_report(left)
    record = right["records"][0]
    record["result_state"] = "failed-controlled" if record["result_state"] != "failed-controlled" else "applied"
    record["stable_reason"] = "backend-failure"
    record["selected_source"] = "artificial-source"
    record["candidate_sources"] = ["artificial-source"]
    record["metrics"]["before"]["length_rms_error"] += metric_tolerance("length_rms_error") * 100.0
    record["snapshot"] = {"changed": True}
    record["policy_evidence"] = {"changed": True}

    diff = compare_baseline_reports(left, right)
    changed_fields = {field["field"] for field in diff["changed_cases"][0]["fields"]}

    assert diff["equivalent"] is False
    assert diff["observational_only"] is True
    assert {
        "result_state",
        "stable_reason",
        "selected_source",
        "candidate_sources",
        "metrics",
        "snapshot",
        "policy_evidence",
    } <= changed_fields


def test_report_diff_ignores_metric_changes_within_tolerance() -> None:
    left = canonicalize_baseline_report(build_baseline_report())
    right = canonicalize_baseline_report(left)
    right["records"][0]["metrics"]["before"]["length_rms_error"] += metric_tolerance("length_rms_error") / 2.0

    diff = compare_baseline_reports(left, right)

    assert diff["equivalent"] is True
    assert diff["changed_cases"] == []


def test_report_diff_detects_added_and_removed_cases() -> None:
    left = canonicalize_baseline_report(build_baseline_report())
    right = canonicalize_baseline_report(left)
    removed = right["records"].pop(0)["case_name"]
    extra = dict(right["records"][0])
    extra["case_name"] = "artificial_extra_case"
    right["records"].append(extra)

    diff = compare_baseline_reports(left, right)

    assert removed in diff["removed_cases"]
    assert "artificial_extra_case" in diff["added_cases"]
    assert diff["observational_only"] is True


def test_baseline_reports_are_test_only() -> None:
    import tests.clean2d_regression.reports as reports

    assert "/tests/clean2d_regression/" in reports.__file__
