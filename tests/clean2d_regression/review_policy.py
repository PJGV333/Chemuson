from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any, Mapping


FIELD_POLICY: dict[str, tuple[str, str]] = {
    "result_state": ("contract", "contract-risk"),
    "stable_reason": ("contract", "needs-review"),
    "selected_source": ("routing", "needs-review"),
    "candidate_sources": ("routing", "needs-review"),
    "metrics": ("geometry-diagnostic", "geometry-risk"),
    "snapshot": ("diagnostic", "needs-review"),
    "policy_evidence": ("complex-policy", "needs-review"),
}


@dataclass(frozen=True)
class Clean2DDiffReviewItem:
    case_name: str
    field: str
    category: str
    severity: str
    known_delicate_related: bool
    observational_only: bool
    message: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify_baseline_diff(
    diff: Mapping[str, Any],
    *,
    report: Mapping[str, Any] | None = None,
) -> tuple[Clean2DDiffReviewItem, ...]:
    known_delicate_cases = _known_delicate_cases(report)
    items = []
    for changed_case in diff.get("changed_cases", ()) or ():
        case_name = str(changed_case.get("case_name", ""))
        known_delicate = case_name in known_delicate_cases
        for changed_field in changed_case.get("fields", ()) or ():
            field = str(changed_field.get("field", ""))
            category, severity = FIELD_POLICY.get(field, ("unknown", "needs-review"))
            items.append(
                Clean2DDiffReviewItem(
                    case_name=case_name,
                    field=field,
                    category=category,
                    severity="known-delicate-change" if known_delicate else severity,
                    known_delicate_related=known_delicate,
                    observational_only=True,
                    message=_message(case_name, field, category, severity, known_delicate),
                )
            )
    json.dumps([item.as_dict() for item in items], allow_nan=False, sort_keys=True)
    return tuple(items)


def summarize_baseline_diff_review(items: tuple[Clean2DDiffReviewItem, ...]) -> dict[str, Any]:
    categories: dict[str, int] = {}
    severities: dict[str, int] = {}
    known_delicate_count = 0
    for item in items:
        categories[item.category] = categories.get(item.category, 0) + 1
        severities[item.severity] = severities.get(item.severity, 0) + 1
        if item.known_delicate_related:
            known_delicate_count += 1
    summary = {
        "equivalent": len(items) == 0,
        "observational_only": True,
        "item_count": len(items),
        "categories": dict(sorted(categories.items())),
        "severities": dict(sorted(severities.items())),
        "known_delicate_related_count": known_delicate_count,
        "items": [item.as_dict() for item in items],
    }
    json.dumps(summary, allow_nan=False, sort_keys=True)
    return summary


def _known_delicate_cases(report: Mapping[str, Any] | None) -> set[str]:
    if report is None:
        return set()
    cases = set()
    for record in report.get("records", ()) or ():
        if "known_delicate" in (record.get("tags", ()) or ()):
            cases.add(str(record.get("case_name", "")))
    return cases


def _message(case_name: str, field: str, category: str, severity: str, known_delicate: bool) -> str:
    prefix = "Known-delicate case changed" if known_delicate else "Baseline diff changed"
    return f"{prefix}: {case_name}.{field} ({category}, {severity})"
