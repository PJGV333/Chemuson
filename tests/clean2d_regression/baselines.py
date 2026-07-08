from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Any, Mapping

from .assertions import execute_case
from .cases import Clean2DRegressionCase
from .metrics import metric_values_equal, validate_before_after_metric_record


SNAPSHOT_EXCLUDED_KEYS = frozenset({"path", "timestamp", "created_at", "absolute_path"})
SNAPSHOT_FLOAT_PRECISION = 9


@dataclass(frozen=True)
class Clean2DBaselineRecord:
    case_name: str
    family: str
    tags: tuple[str, ...]
    mode: str
    target: str
    result_state: str
    stable_reason: str | None
    selected_source: str | None
    candidate_sources: tuple[str, ...]
    metrics: dict[str, Any]
    snapshot: dict[str, Any] | None
    policy_evidence: dict[str, Any] | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "case_name": self.case_name,
            "family": self.family,
            "tags": list(self.tags),
            "mode": self.mode,
            "target": self.target,
            "result_state": self.result_state,
            "stable_reason": self.stable_reason,
            "selected_source": self.selected_source,
            "candidate_sources": list(self.candidate_sources),
            "metrics": self.metrics,
            "snapshot": self.snapshot,
            "policy_evidence": self.policy_evidence,
        }


def build_baseline_record(case: Clean2DRegressionCase) -> Clean2DBaselineRecord:
    execution = execute_case(case)
    snapshot = execution["engine_debug_snapshot"]
    candidates = snapshot.get("candidates", ()) if isinstance(snapshot, Mapping) else ()
    selected_source = _selected_source(candidates)

    record = Clean2DBaselineRecord(
        case_name=str(execution["case"]["name"]),
        family=str(execution["case"]["family"]),
        tags=tuple(str(tag) for tag in execution["case"]["tags"]),
        mode=str(execution["case"]["mode"]),
        target=str(execution["case"]["target"]),
        result_state=str(execution["result"]["state"]),
        stable_reason=execution["result"].get("reason"),
        selected_source=selected_source,
        candidate_sources=tuple(str(source) for source in execution["result"]["candidate_sources"]),
        metrics=validate_before_after_metric_record(execution["metrics"]),
        snapshot=canonicalize_snapshot(snapshot),
        policy_evidence=_policy_evidence(candidates),
    )
    json.dumps(record.as_dict(), allow_nan=False, sort_keys=True)
    return record


def canonicalize_baseline_record(record: Clean2DBaselineRecord | Mapping[str, Any]) -> dict[str, Any]:
    data = record.as_dict() if isinstance(record, Clean2DBaselineRecord) else dict(record)
    canonical = _canonicalize(data)
    json.dumps(canonical, allow_nan=False, sort_keys=True)
    return canonical


def canonicalize_snapshot(snapshot: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if snapshot is None:
        return None
    return _canonicalize(snapshot, exclude_keys=SNAPSHOT_EXCLUDED_KEYS)


def baseline_records_equivalent(
    left: Clean2DBaselineRecord | Mapping[str, Any],
    right: Clean2DBaselineRecord | Mapping[str, Any],
) -> bool:
    left_data = canonicalize_baseline_record(left)
    right_data = canonicalize_baseline_record(right)
    return _equivalent(left_data, right_data, metric_path=())


def _selected_source(candidates: Any) -> str | None:
    if not isinstance(candidates, list):
        return None
    for candidate in candidates:
        if isinstance(candidate, Mapping) and candidate.get("selected") is True:
            source = candidate.get("source")
            return None if source is None else str(source)
    return None


def _policy_evidence(candidates: Any) -> dict[str, Any] | None:
    evidence: dict[str, Any] = {}
    if not isinstance(candidates, list):
        return None
    for candidate in candidates:
        if not isinstance(candidate, Mapping):
            continue
        diagnostic = candidate.get("quality_diagnostic")
        if not isinstance(diagnostic, Mapping):
            continue
        internal = diagnostic.get("internal")
        if not isinstance(internal, Mapping):
            continue
        for key, value in internal.items():
            if str(key).startswith("complex_"):
                evidence[str(key)] = value
    return _canonicalize(evidence) if evidence else None


def _canonicalize(value: Any, *, exclude_keys: frozenset[str] = frozenset()) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise AssertionError("baseline records must not contain NaN or infinity")
        return round(float(value), SNAPSHOT_FLOAT_PRECISION)
    if isinstance(value, Mapping):
        return {
            str(key): _canonicalize(item, exclude_keys=exclude_keys)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            if str(key) not in exclude_keys
        }
    if isinstance(value, (set, frozenset)):
        return [_canonicalize(item, exclude_keys=exclude_keys) for item in sorted(value, key=repr)]
    if isinstance(value, (list, tuple)):
        return [_canonicalize(item, exclude_keys=exclude_keys) for item in value]
    raise AssertionError(f"unstable baseline value: {type(value).__name__}")


def _equivalent(left: Any, right: Any, *, metric_path: tuple[str, ...]) -> bool:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        if set(left) != set(right):
            return False
        return all(
            _equivalent(left[key], right[key], metric_path=(*metric_path, str(key)))
            for key in sorted(left)
        )
    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            return False
        return all(_equivalent(a, b, metric_path=metric_path) for a, b in zip(left, right))
    metric_name = metric_path[-1] if metric_path else ""
    if _inside_metrics(metric_path) and isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return metric_values_equal(metric_name, left, right)
    return left == right


def _inside_metrics(path: tuple[str, ...]) -> bool:
    return len(path) >= 3 and path[-3] == "metrics" and path[-2] in {"before", "after"}
