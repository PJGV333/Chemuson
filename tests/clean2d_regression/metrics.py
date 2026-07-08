from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Any, Mapping


VALUE_TYPES = frozenset({"string", "integer", "number", "nullable-number"})
POLARITIES = frozenset({"lower-is-better", "higher-is-better", "categorical", "informational"})
SCOPES = frozenset({"selected-subgraph", "final-result"})

NORMALIZED_FLOAT_TOLERANCE = 1e-9
PIXEL_FLOAT_TOLERANCE = 1e-6
DEGREE_FLOAT_TOLERANCE = 1e-6


@dataclass(frozen=True)
class Clean2DMetricDefinition:
    name: str
    value_type: str
    unit: str | None
    polarity: str
    scope: str
    required: bool
    diagnostic_only: bool
    optional_when: str = ""


CLEAN2D_GEOMETRY_METRICS: dict[str, Clean2DMetricDefinition] = {
    "quality_class": Clean2DMetricDefinition(
        name="quality_class",
        value_type="string",
        unit=None,
        polarity="categorical",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "reason": Clean2DMetricDefinition(
        name="reason",
        value_type="string",
        unit=None,
        polarity="informational",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
        optional_when="No quality reason applies.",
    ),
    "crossings": Clean2DMetricDefinition(
        name="crossings",
        value_type="integer",
        unit="count",
        polarity="lower-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "min_nonbonded_distance": Clean2DMetricDefinition(
        name="min_nonbonded_distance",
        value_type="nullable-number",
        unit="canvas-px",
        polarity="higher-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
        optional_when="No finite non-bonded atom pair exists in the evaluated scope.",
    ),
    "min_ring_degeneracy": Clean2DMetricDefinition(
        name="min_ring_degeneracy",
        value_type="nullable-number",
        unit="dimensionless",
        polarity="higher-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
        optional_when="No ring exists in the evaluated scope.",
    ),
    "length_rms_error": Clean2DMetricDefinition(
        name="length_rms_error",
        value_type="number",
        unit="normalized-to-target-bond-length",
        polarity="lower-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "length_max_error": Clean2DMetricDefinition(
        name="length_max_error",
        value_type="number",
        unit="normalized-to-target-bond-length",
        polarity="lower-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "angle_rms_deviation": Clean2DMetricDefinition(
        name="angle_rms_deviation",
        value_type="number",
        unit="degrees",
        polarity="lower-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "angle_max_deviation": Clean2DMetricDefinition(
        name="angle_max_deviation",
        value_type="number",
        unit="degrees",
        polarity="lower-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
    "visual_score": Clean2DMetricDefinition(
        name="visual_score",
        value_type="number",
        unit="normalized-score",
        polarity="higher-is-better",
        scope="selected-subgraph",
        required=True,
        diagnostic_only=True,
    ),
}


def validate_metric_definitions() -> None:
    for name, definition in CLEAN2D_GEOMETRY_METRICS.items():
        if definition.name != name:
            raise AssertionError(f"metric definition key/name mismatch: {name}")
        if definition.value_type not in VALUE_TYPES:
            raise AssertionError(f"invalid value_type for {name}: {definition.value_type}")
        if definition.polarity not in POLARITIES:
            raise AssertionError(f"invalid polarity for {name}: {definition.polarity}")
        if definition.scope not in SCOPES:
            raise AssertionError(f"invalid scope for {name}: {definition.scope}")
        if not definition.diagnostic_only:
            raise AssertionError(f"metric must remain diagnostic-only in this change: {name}")


def validate_metric_record(record: Mapping[str, Any]) -> dict[str, Any]:
    normalized = {str(key): _json_metric_value(key, value) for key, value in record.items()}
    missing = [
        name
        for name, definition in CLEAN2D_GEOMETRY_METRICS.items()
        if definition.required and name not in normalized
    ]
    if missing:
        raise AssertionError(f"missing required metrics: {missing}")
    extra = sorted(set(normalized) - set(CLEAN2D_GEOMETRY_METRICS))
    if extra:
        raise AssertionError(f"undefined metrics emitted: {extra}")
    for name, value in normalized.items():
        _validate_metric_value(CLEAN2D_GEOMETRY_METRICS[name], value)
    json.dumps(normalized, allow_nan=False, sort_keys=True)
    return normalized


def validate_before_after_metric_record(record: Mapping[str, Any]) -> dict[str, Any]:
    before = validate_metric_record(_mapping_value(record, "before"))
    after_value = record.get("after")
    after = None if after_value is None else validate_metric_record(_mapping_value(record, "after"))
    result = {"before": before, "after": after}
    json.dumps(result, allow_nan=False, sort_keys=True)
    return result


def metric_values_equal(name: str, left: Any, right: Any) -> bool:
    definition = CLEAN2D_GEOMETRY_METRICS[name]
    if definition.value_type not in {"number", "nullable-number"}:
        return left == right
    if left is None or right is None:
        return left is None and right is None
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=metric_tolerance(name))


def metric_tolerance(name: str) -> float:
    unit = CLEAN2D_GEOMETRY_METRICS[name].unit
    if unit == "canvas-px":
        return PIXEL_FLOAT_TOLERANCE
    if unit == "degrees":
        return DEGREE_FLOAT_TOLERANCE
    return NORMALIZED_FLOAT_TOLERANCE


def _mapping_value(record: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = record.get(key)
    if not isinstance(value, Mapping):
        raise AssertionError(f"metric record field must be a mapping: {key}")
    return value


def _json_metric_value(name: str, value: Any) -> Any:
    if isinstance(value, bool):
        raise AssertionError(f"boolean is not a valid metric value for {name}")
    if isinstance(value, float):
        if math.isfinite(value):
            return float(value)
        return None
    if isinstance(value, int):
        return int(value)
    if value is None or isinstance(value, str):
        return value
    raise AssertionError(f"unstable metric value for {name}: {type(value).__name__}")


def _validate_metric_value(definition: Clean2DMetricDefinition, value: Any) -> None:
    if value is None:
        if definition.value_type != "nullable-number" and not definition.optional_when:
            raise AssertionError(f"metric {definition.name} does not allow None")
        return
    if definition.value_type == "string":
        if not isinstance(value, str):
            raise AssertionError(f"metric {definition.name} must be string")
        return
    if definition.value_type == "integer":
        if not isinstance(value, int) or isinstance(value, bool):
            raise AssertionError(f"metric {definition.name} must be integer")
        if value < 0:
            raise AssertionError(f"metric {definition.name} must be non-negative")
        return
    if definition.value_type in {"number", "nullable-number"}:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise AssertionError(f"metric {definition.name} must be numeric")
        if not math.isfinite(float(value)):
            raise AssertionError(f"metric {definition.name} must be finite or None")
        return
    raise AssertionError(f"unknown value type for {definition.name}: {definition.value_type}")
