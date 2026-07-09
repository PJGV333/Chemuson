from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping


BLOCKING_CANDIDATE = "blocking-candidate"
REVIEW_REQUIRED = "review-required"
INFORMATIONAL = "informational"
ALLOWED_KNOWN_DELICATE = "allowed-known-delicate"
FUTURE_GATE_CANDIDATE = "future-gate-candidate"
CLEAN = "clean"

_REVIEW_REQUIRED_SEVERITIES = {"needs-review", "geometry-risk"}


@dataclass(frozen=True)
class Clean2DDiffBlockingDecision:
    decision: str
    blocking_candidate_count: int
    review_required_count: int
    informational_count: int
    allowed_known_delicate_count: int
    future_gate_candidate_count: int
    reasons: tuple[dict[str, Any], ...]
    observational_only: bool

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def evaluate_diff_blocking_policy(review_summary: Mapping[str, Any]) -> Clean2DDiffBlockingDecision:
    reasons = tuple(_reason_for_item(item) for item in _review_items(review_summary))
    blocking_candidate_count = sum(
        1 for reason in reasons if reason["policy_classification"] == BLOCKING_CANDIDATE
    )
    review_required_count = sum(1 for reason in reasons if reason["policy_classification"] == REVIEW_REQUIRED)
    informational_count = sum(1 for reason in reasons if reason["policy_classification"] == INFORMATIONAL)
    allowed_known_delicate_count = sum(
        1 for reason in reasons if reason["policy_classification"] == ALLOWED_KNOWN_DELICATE
    )
    future_gate_candidate_count = sum(
        1 for reason in reasons if reason["policy_classification"] == FUTURE_GATE_CANDIDATE
    )

    if blocking_candidate_count:
        decision = BLOCKING_CANDIDATE
    elif review_required_count:
        decision = REVIEW_REQUIRED
    elif allowed_known_delicate_count and allowed_known_delicate_count == len(reasons):
        decision = ALLOWED_KNOWN_DELICATE
    elif future_gate_candidate_count:
        decision = FUTURE_GATE_CANDIDATE
    elif informational_count:
        decision = INFORMATIONAL
    else:
        decision = CLEAN

    result = Clean2DDiffBlockingDecision(
        decision=decision,
        blocking_candidate_count=blocking_candidate_count,
        review_required_count=review_required_count,
        informational_count=informational_count,
        allowed_known_delicate_count=allowed_known_delicate_count,
        future_gate_candidate_count=future_gate_candidate_count,
        reasons=reasons,
        observational_only=True,
    )
    json.dumps(result.as_dict(), allow_nan=False, sort_keys=True)
    return result


def _review_items(review_summary: Mapping[str, Any]) -> Iterable[Mapping[str, Any]]:
    items = review_summary.get("items", ()) or ()
    for item in items:
        if isinstance(item, Mapping):
            yield item


def _reason_for_item(item: Mapping[str, Any]) -> dict[str, Any]:
    severity = str(item.get("severity") or item.get("decision") or item.get("reason") or INFORMATIONAL)
    field = str(item.get("field", ""))
    known_delicate = bool(item.get("known_delicate_related") or item.get("known_delicate"))
    classification = _classification_for(severity, known_delicate)
    return {
        "case_name": str(item.get("case_name", "")),
        "category": str(item.get("category", "")),
        "field": field,
        "input_severity": severity,
        "known_delicate_related": known_delicate or severity == "known-delicate-change",
        "message": str(item.get("message", "")),
        "observational_only": True,
        "policy_classification": classification,
    }


def _classification_for(severity: str, known_delicate: bool) -> str:
    if severity == "known-delicate-change" or known_delicate:
        return ALLOWED_KNOWN_DELICATE
    if severity == "contract-risk":
        return BLOCKING_CANDIDATE
    if severity in _REVIEW_REQUIRED_SEVERITIES:
        return REVIEW_REQUIRED
    if severity == FUTURE_GATE_CANDIDATE:
        return FUTURE_GATE_CANDIDATE
    return INFORMATIONAL
