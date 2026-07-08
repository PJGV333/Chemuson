from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Mapping

from .baselines import (
    Clean2DBaselineRecord,
    baseline_records_equivalent,
    build_baseline_record,
    canonicalize_baseline_record,
)
from .cases import Clean2DRegressionCase, get_regression_cases


CLEAN2D_BASELINE_REPORT_SCHEMA = "chemuson.clean2d.baseline-report"
CLEAN2D_BASELINE_REPORT_VERSION = 1
OBSERVABLE_DIFF_FIELDS = (
    "result_state",
    "stable_reason",
    "selected_source",
    "candidate_sources",
    "metrics",
    "snapshot",
    "policy_evidence",
)


@dataclass(frozen=True)
class Clean2DBaselineReport:
    schema: str
    version: int
    metadata: dict[str, Any]
    records: tuple[dict[str, Any], ...]
    summary: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "version": self.version,
            "metadata": self.metadata,
            "summary": self.summary,
            "records": list(self.records),
        }


def build_baseline_report(
    cases: tuple[Clean2DRegressionCase, ...] | None = None,
    *,
    metadata: Mapping[str, Any] | None = None,
) -> Clean2DBaselineReport:
    selected_cases = tuple(cases or get_regression_cases())
    records = tuple(
        sorted(
            (canonicalize_baseline_record(build_baseline_record(case)) for case in selected_cases),
            key=lambda record: record["case_name"],
        )
    )
    report = Clean2DBaselineReport(
        schema=CLEAN2D_BASELINE_REPORT_SCHEMA,
        version=CLEAN2D_BASELINE_REPORT_VERSION,
        metadata=_canonicalize(metadata or {}),
        records=records,
        summary=_summary(records),
    )
    json.dumps(report.as_dict(), allow_nan=False, sort_keys=True)
    return report


def canonicalize_baseline_report(report: Clean2DBaselineReport | Mapping[str, Any]) -> dict[str, Any]:
    data = report.as_dict() if isinstance(report, Clean2DBaselineReport) else dict(report)
    records = tuple(sorted((_canonicalize(record) for record in data.get("records", ())), key=lambda row: row["case_name"]))
    canonical = {
        "schema": str(data.get("schema", CLEAN2D_BASELINE_REPORT_SCHEMA)),
        "version": int(data.get("version", CLEAN2D_BASELINE_REPORT_VERSION)),
        "metadata": _canonicalize(data.get("metadata", {})),
        "summary": _canonicalize(data.get("summary", _summary(records))),
        "records": list(records),
    }
    json.dumps(canonical, allow_nan=False, sort_keys=True)
    return canonical


def baseline_report_to_json(report: Clean2DBaselineReport | Mapping[str, Any]) -> str:
    return json.dumps(canonicalize_baseline_report(report), allow_nan=False, indent=2, sort_keys=True) + "\n"


def compare_baseline_reports(
    left: Clean2DBaselineReport | Mapping[str, Any],
    right: Clean2DBaselineReport | Mapping[str, Any],
) -> dict[str, Any]:
    left_report = canonicalize_baseline_report(left)
    right_report = canonicalize_baseline_report(right)
    left_records = {record["case_name"]: record for record in left_report["records"]}
    right_records = {record["case_name"]: record for record in right_report["records"]}
    added = sorted(set(right_records) - set(left_records))
    removed = sorted(set(left_records) - set(right_records))
    changed = []
    for case_name in sorted(set(left_records) & set(right_records)):
        left_record = left_records[case_name]
        right_record = right_records[case_name]
        if baseline_records_equivalent(left_record, right_record):
            continue
        fields = []
        for field in OBSERVABLE_DIFF_FIELDS:
            probe_left = dict(left_record)
            probe_right = dict(right_record)
            probe_left[field] = right_record.get(field)
            if not baseline_records_equivalent(probe_left, right_record):
                fields.append(
                    {
                        "field": field,
                        "before": left_record.get(field),
                        "after": right_record.get(field),
                    }
                )
        changed.append({"case_name": case_name, "fields": fields})
    diff = {
        "equivalent": not added and not removed and not changed,
        "observational_only": True,
        "added_cases": added,
        "removed_cases": removed,
        "changed_cases": changed,
    }
    json.dumps(diff, allow_nan=False, sort_keys=True)
    return diff


def _summary(records: tuple[Mapping[str, Any], ...]) -> dict[str, Any]:
    states: dict[str, int] = {}
    tags: dict[str, int] = {}
    families: dict[str, int] = {}
    for record in records:
        states[str(record["result_state"])] = states.get(str(record["result_state"]), 0) + 1
        families[str(record["family"])] = families.get(str(record["family"]), 0) + 1
        for tag in record.get("tags", ()):
            tags[str(tag)] = tags.get(str(tag), 0) + 1
    return {
        "case_count": len(records),
        "result_states": dict(sorted(states.items())),
        "families": dict(sorted(families.items())),
        "tags": dict(sorted(tags.items())),
    }


def _canonicalize(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _canonicalize(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, (set, frozenset)):
        return [_canonicalize(item) for item in sorted(value, key=repr)]
    if isinstance(value, (list, tuple)):
        return [_canonicalize(item) for item in value]
    if isinstance(value, Clean2DBaselineRecord):
        return canonicalize_baseline_record(value)
    raise AssertionError(f"unstable report value: {type(value).__name__}")
