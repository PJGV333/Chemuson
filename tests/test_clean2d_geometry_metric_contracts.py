from __future__ import annotations

import json

from tests.clean2d_regression.assertions import execute_case
from tests.clean2d_regression.cases import get_regression_cases
from tests.clean2d_regression.metrics import (
    CLEAN2D_GEOMETRY_METRICS,
    metric_tolerance,
    metric_values_equal,
    validate_before_after_metric_record,
    validate_metric_definitions,
    validate_metric_record,
)


EXPECTED_METRICS = {
    "quality_class",
    "reason",
    "crossings",
    "min_nonbonded_distance",
    "min_ring_degeneracy",
    "length_rms_error",
    "length_max_error",
    "angle_rms_deviation",
    "angle_max_deviation",
    "visual_score",
}


def test_geometry_metric_registry_defines_all_corpus_metrics() -> None:
    validate_metric_definitions()

    assert set(CLEAN2D_GEOMETRY_METRICS) == EXPECTED_METRICS
    for name, definition in CLEAN2D_GEOMETRY_METRICS.items():
        assert definition.name == name
        assert definition.value_type
        assert definition.polarity
        assert definition.scope
        assert definition.diagnostic_only is True


def test_metric_definitions_capture_units_and_polarity() -> None:
    assert CLEAN2D_GEOMETRY_METRICS["quality_class"].polarity == "categorical"
    assert CLEAN2D_GEOMETRY_METRICS["reason"].polarity == "informational"
    assert CLEAN2D_GEOMETRY_METRICS["crossings"].unit == "count"
    assert CLEAN2D_GEOMETRY_METRICS["crossings"].polarity == "lower-is-better"
    assert CLEAN2D_GEOMETRY_METRICS["min_nonbonded_distance"].unit == "canvas-px"
    assert CLEAN2D_GEOMETRY_METRICS["min_nonbonded_distance"].polarity == "higher-is-better"
    assert CLEAN2D_GEOMETRY_METRICS["min_ring_degeneracy"].unit == "dimensionless"
    assert CLEAN2D_GEOMETRY_METRICS["visual_score"].polarity == "higher-is-better"
    assert CLEAN2D_GEOMETRY_METRICS["angle_rms_deviation"].unit == "degrees"


def test_metric_records_serialize_as_json_stable_primitives() -> None:
    record = {
        "quality_class": "good",
        "reason": "",
        "crossings": 0,
        "min_nonbonded_distance": None,
        "min_ring_degeneracy": None,
        "length_rms_error": 0.0,
        "length_max_error": 0.0,
        "angle_rms_deviation": 0.0,
        "angle_max_deviation": 0.0,
        "visual_score": 1.0,
    }

    normalized = validate_metric_record(record)

    assert normalized == record
    json.dumps(normalized, allow_nan=False, sort_keys=True)


def test_float_metric_comparison_uses_explicit_tolerances() -> None:
    assert metric_tolerance("length_rms_error") == 1e-9
    assert metric_tolerance("min_nonbonded_distance") == 1e-6
    assert metric_tolerance("angle_max_deviation") == 1e-6
    assert metric_values_equal("length_rms_error", 1.0, 1.0 + 5e-10)
    assert not metric_values_equal("length_rms_error", 1.0, 1.0 + 5e-8)


def test_optional_metrics_may_be_none_only_when_declared() -> None:
    record = {
        "quality_class": "good",
        "reason": None,
        "crossings": 0,
        "min_nonbonded_distance": None,
        "min_ring_degeneracy": None,
        "length_rms_error": 0.0,
        "length_max_error": 0.0,
        "angle_rms_deviation": 0.0,
        "angle_max_deviation": 0.0,
        "visual_score": 1.0,
    }

    normalized = validate_metric_record(record)

    assert normalized["reason"] is None
    assert normalized["min_nonbonded_distance"] is None
    assert normalized["min_ring_degeneracy"] is None


def test_expanded_corpus_emits_metric_records_compatible_with_contract() -> None:
    for case in get_regression_cases():
        snapshot = execute_case(case)
        metrics = validate_before_after_metric_record(snapshot["metrics"])

        assert set(metrics["before"]) == EXPECTED_METRICS
        if metrics["after"] is not None:
            assert set(metrics["after"]) == EXPECTED_METRICS


def test_metric_contract_remains_diagnostic_only_for_known_delicate_cases() -> None:
    delicate = [case for case in get_regression_cases() if case.has_tag("known_delicate")]
    assert delicate
    for case in delicate:
        snapshot = execute_case(case)
        validate_before_after_metric_record(snapshot["metrics"])
        assert snapshot["result"]["state"] in case.expected_states


def test_metric_contract_does_not_recompute_result_state_from_values() -> None:
    case = get_regression_cases()[0]
    snapshot = execute_case(case)
    before_state = snapshot["result"]["state"]

    validate_before_after_metric_record(snapshot["metrics"])

    assert snapshot["result"]["state"] == before_state
