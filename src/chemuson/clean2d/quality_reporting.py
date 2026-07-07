from __future__ import annotations

"""Stable reporting helpers for Clean 2D quality diagnostics.

These helpers are intentionally reporting-only: they validate and normalize
diagnostic output without changing layout, ranking, or candidate selection.
"""

from collections.abc import Mapping
from math import isfinite
from typing import Any, Literal, TypedDict, cast


Clean2DQualityState = Literal[
    "applied",
    "rejected",
    "preserve-only",
    "no-op",
    "failed-controlled",
]

Clean2DQualityReason = Literal[
    "invalid-coordinates",
    "invariant-violation",
    "stereo-risk",
    "boundary-bond-risk",
    "new-crossing-risk",
    "collision-risk",
    "collapsed-ring-risk",
    "excessive-displacement",
    "worse-quality",
    "backend-failure",
]


class Clean2DQualityDiagnostic(TypedDict, total=False):
    state: Clean2DQualityState
    reason: Clean2DQualityReason
    score: float
    source: str
    internal: dict[str, Any]


class Clean2DQualityContractError(ValueError):
    """Raised when a Clean 2D quality diagnostic violates the stable contract."""


STABLE_CLEAN2D_QUALITY_STATES: frozenset[str] = frozenset(
    {"applied", "rejected", "preserve-only", "no-op", "failed-controlled"}
)

STABLE_CLEAN2D_QUALITY_REASONS: frozenset[str] = frozenset(
    {
        "invalid-coordinates",
        "invariant-violation",
        "stereo-risk",
        "boundary-bond-risk",
        "new-crossing-risk",
        "collision-risk",
        "collapsed-ring-risk",
        "excessive-displacement",
        "worse-quality",
        "backend-failure",
    }
)


def normalize_clean2d_reporting_score(
    raw_score: float | int | None,
    *,
    lower_is_better: bool = True,
    scale: float = 100.0,
) -> float:
    """Normalize a diagnostic score to 0-1 for reporting only.

    This function does not feed back into ranking or selection. For
    lower-is-better raw scores, 0 maps to 1.0 and ``scale`` maps to 0.0.
    """
    try:
        raw = float(raw_score if raw_score is not None else 0.0)
    except (TypeError, ValueError):
        raw = scale
    if not isfinite(raw):
        raw = scale
    safe_scale = max(float(scale or 1.0), 1e-9)
    if lower_is_better:
        normalized = 1.0 - (raw / safe_scale)
    else:
        normalized = raw / safe_scale
    return max(0.0, min(1.0, normalized))


def clean2d_quality_diagnostic(
    *,
    state: str,
    source: str,
    score: float | int,
    reason: str | None = None,
    internal: Mapping[str, Any] | None = None,
) -> Clean2DQualityDiagnostic:
    diagnostic: Clean2DQualityDiagnostic = {
        "state": cast(Clean2DQualityState, state),
        "score": float(score),
        "source": str(source or "").strip(),
    }
    if reason is not None and str(reason).strip():
        diagnostic["reason"] = cast(Clean2DQualityReason, str(reason).strip())
    if internal is not None:
        diagnostic["internal"] = dict(internal)
    validate_clean2d_quality_diagnostic(diagnostic)
    return diagnostic


def validate_clean2d_quality_diagnostic(diagnostic: Mapping[str, Any]) -> None:
    state = str(diagnostic.get("state", "")).strip()
    if state not in STABLE_CLEAN2D_QUALITY_STATES:
        raise Clean2DQualityContractError(f"invalid Clean 2D diagnostic state: {state!r}")

    source = str(diagnostic.get("source", "")).strip()
    if not source:
        raise Clean2DQualityContractError("Clean 2D diagnostic source is required")

    try:
        score = float(diagnostic["score"])
    except (KeyError, TypeError, ValueError) as exc:
        raise Clean2DQualityContractError("Clean 2D diagnostic score must be numeric") from exc
    if not isfinite(score) or not 0.0 <= score <= 1.0:
        raise Clean2DQualityContractError("Clean 2D diagnostic score must be in [0, 1]")

    reason = diagnostic.get("reason")
    has_reason = reason is not None and str(reason).strip() != ""
    if state in {"rejected", "failed-controlled"} and not has_reason:
        raise Clean2DQualityContractError(f"Clean 2D diagnostic reason is required for {state}")
    if state == "applied" and has_reason:
        raise Clean2DQualityContractError("Clean 2D applied diagnostics must omit contractual reason")
    if has_reason and str(reason) not in STABLE_CLEAN2D_QUALITY_REASONS:
        raise Clean2DQualityContractError(f"invalid Clean 2D diagnostic reason: {reason!r}")

    internal = diagnostic.get("internal")
    if internal is not None and not isinstance(internal, dict):
        raise Clean2DQualityContractError("Clean 2D diagnostic internal metadata must be a dict")


def clean2d_safety_quality_diagnostic(report: Any) -> Clean2DQualityDiagnostic:
    """Adapt an existing safety report surface to the stable diagnostic contract."""
    state = "applied" if bool(getattr(report, "passed", False)) else "rejected"
    raw_reason = str(getattr(report, "rejection_reason", "") or "")
    reason = _stable_clean2d_reason(raw_reason) if state == "rejected" else None
    score = _safety_reporting_score(report)
    return clean2d_quality_diagnostic(
        state=state,
        source="safety",
        score=score,
        reason=reason,
        internal={
            "rejection_reason": raw_reason,
            "missing_coords": tuple(getattr(report, "missing_coords", ()) or ()),
            "nan_coords": tuple(getattr(report, "nan_coords", ()) or ()),
            "new_crossings": int(getattr(report, "new_crossings", 0) or 0),
            "min_nonbonded_after": float(getattr(report, "min_nonbonded_after", 0.0) or 0.0),
            "ring_degeneracy_after": float(getattr(report, "ring_degeneracy_after", 0.0) or 0.0),
            "max_displacement": float(getattr(report, "max_displacement", 0.0) or 0.0),
            "target_bond_length": float(getattr(report, "target_bond_length", 0.0) or 0.0),
        },
    )


def clean2d_candidate_quality_diagnostic(candidate: Any) -> Clean2DQualityDiagnostic:
    state = str(getattr(candidate, "outcome_state", "rejected"))
    reason = None
    if state == "rejected":
        reason = str(getattr(candidate, "stable_rejection_reason", "backend-failure"))
    raw_score = getattr(candidate, "score", 0.0)
    return clean2d_quality_diagnostic(
        state=state,
        source=str(getattr(candidate, "source", "candidate") or "candidate"),
        score=normalize_clean2d_reporting_score(raw_score, lower_is_better=True, scale=100.0),
        reason=reason,
        internal={
            "raw_score": raw_score,
            "novelty": getattr(candidate, "novelty", 0.0),
            "rejection_reason": getattr(candidate, "rejection_reason", ""),
            "geometry_hash": getattr(candidate, "geometry_hash", ""),
            "metadata": dict(getattr(candidate, "metadata", {}) or {}),
        },
    )


def clean2d_result_quality_diagnostic(result: Any) -> Clean2DQualityDiagnostic:
    selected = getattr(result, "selected", None)
    rejected = tuple(getattr(result, "rejected", ()) or ())
    candidates = tuple(getattr(result, "candidates", ()) or ())
    source = "engine"
    raw_score: float | int | None = 100.0
    if selected is not None:
        source = str(getattr(selected, "source", "engine") or "engine")
        raw_score = getattr(selected, "score", 0.0)
    elif rejected:
        source = str(getattr(rejected[0], "source", "engine") or "engine")
        raw_score = getattr(rejected[0], "score", 100.0)

    state = str(getattr(result, "result_state", "failed-controlled"))
    reason = None
    if state in {"rejected", "failed-controlled"}:
        reason = str(getattr(result, "stable_reason", "backend-failure") or "backend-failure")
    return clean2d_quality_diagnostic(
        state=state,
        source=source,
        score=normalize_clean2d_reporting_score(raw_score, lower_is_better=True, scale=100.0),
        reason=reason,
        internal={
            "candidate_sources": tuple(getattr(result, "candidate_sources", ())),
            "candidate_count": len(candidates),
            "rejected_count": len(rejected),
            "selected_score": raw_score if selected is not None else None,
            "message": getattr(result, "message", ""),
        },
    )


def _safety_reporting_score(report: Any) -> float:
    if bool(getattr(report, "passed", False)):
        return 1.0
    target = max(1.0, float(getattr(report, "target_bond_length", 1.0) or 1.0))
    penalties = [0.0]
    penalties.append(float(getattr(report, "new_crossings", 0) or 0) * 50.0)
    min_nonbonded = float(getattr(report, "min_nonbonded_after", target) or target)
    if min_nonbonded < target * 0.25:
        penalties.append((target * 0.25 - min_nonbonded) / target * 100.0)
    ring_degeneracy = float(getattr(report, "ring_degeneracy_after", 1.0) or 0.0)
    if bool(getattr(report, "is_cyclic", False)) and ring_degeneracy < 0.05:
        penalties.append((0.05 - ring_degeneracy) * 100.0)
    max_displacement = float(getattr(report, "max_displacement", 0.0) or 0.0)
    penalties.append(max(0.0, max_displacement / target - 1.0) * 10.0)
    if getattr(report, "missing_coords", None) or getattr(report, "nan_coords", None):
        penalties.append(100.0)
    return normalize_clean2d_reporting_score(max(penalties), lower_is_better=True, scale=100.0)


def _stable_clean2d_reason(reason: str | None) -> str:
    from chemuson.clean2d.engine import stable_clean2d_rejection_reason

    return stable_clean2d_rejection_reason(reason)
