from __future__ import annotations

import json

import pytest

from tests.clean2d_regression.baselines import (
    SNAPSHOT_EXCLUDED_KEYS,
    baseline_records_equivalent,
    build_baseline_record,
    canonicalize_baseline_record,
    canonicalize_snapshot,
)
from tests.clean2d_regression.cases import get_regression_cases
from tests.clean2d_regression.metrics import CLEAN2D_GEOMETRY_METRICS, metric_tolerance


def test_corpus_case_order_is_stable() -> None:
    names = [case.name for case in get_regression_cases()]

    assert names == [case.name for case in get_regression_cases()]
    assert names == list(dict.fromkeys(names))


@pytest.mark.parametrize("case", get_regression_cases(), ids=lambda case: case.name)
def test_baseline_record_is_json_serializable(case) -> None:
    record = build_baseline_record(case)
    data = canonicalize_baseline_record(record)

    json.dumps(data, allow_nan=False, sort_keys=True)
    assert data["case_name"] == case.name
    assert data["family"] == case.family
    assert data["tags"] == list(case.tags)
    assert data["mode"] == case.mode.value
    assert data["target"] == case.target_label
    assert data["result_state"] in case.expected_states
    assert "stable_reason" in data
    assert set(data["metrics"]["before"]) == set(CLEAN2D_GEOMETRY_METRICS)


@pytest.mark.parametrize("case", get_regression_cases(), ids=lambda case: case.name)
def test_repeated_case_execution_produces_equivalent_baseline(case) -> None:
    first = build_baseline_record(case)
    second = build_baseline_record(case)

    assert baseline_records_equivalent(first, second), case.name


def test_baseline_comparison_uses_metric_tolerances() -> None:
    case = get_regression_cases()[0]
    first = canonicalize_baseline_record(build_baseline_record(case))
    second = canonicalize_baseline_record(first)
    second["metrics"]["before"]["length_rms_error"] += metric_tolerance("length_rms_error") / 2.0

    assert baseline_records_equivalent(first, second)

    second["metrics"]["before"]["length_rms_error"] += metric_tolerance("length_rms_error") * 10.0
    assert not baseline_records_equivalent(first, second)


def test_metrics_remain_diagnostic_only_and_result_state_is_not_recomputed() -> None:
    case = next(item for item in get_regression_cases() if item.has_tag("known_delicate"))
    record = canonicalize_baseline_record(build_baseline_record(case))
    result_state = record["result_state"]

    record["metrics"]["before"]["visual_score"] = -999.0

    assert record["result_state"] == result_state
    assert CLEAN2D_GEOMETRY_METRICS["visual_score"].diagnostic_only is True


def test_snapshot_canonicalization_excludes_documented_ephemeral_fields() -> None:
    snapshot = {
        "schema": "x",
        "timestamp": "now",
        "path": "/tmp/example.json",
        "nested": {"absolute_path": "/tmp/other", "keep": 1.0},
    }

    canonical = canonicalize_snapshot(snapshot)

    assert canonical == {"nested": {"keep": 1.0}, "schema": "x"}
    assert {"timestamp", "path", "absolute_path"} <= SNAPSHOT_EXCLUDED_KEYS


@pytest.mark.parametrize("case", get_regression_cases(), ids=lambda case: case.name)
def test_candidate_sources_are_stable_between_repeated_runs(case) -> None:
    first = build_baseline_record(case)
    second = build_baseline_record(case)

    assert first.candidate_sources == second.candidate_sources
    assert first.selected_source == second.selected_source
    assert first.result_state == second.result_state
    assert first.stable_reason == second.stable_reason


def test_baseline_harness_is_test_only() -> None:
    import tests.clean2d_regression.baselines as baselines

    assert "/tests/clean2d_regression/" in baselines.__file__
