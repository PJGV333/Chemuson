"""Opt-in diagnostic evaluation of validated molecular-assistant proposals.

This command is a campaign tool, not an application service. Importing it never
contacts a provider; network I/O occurs only from ``main`` after explicit CLI
configuration.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
from collections.abc import Callable, Sequence
from typing import Any

from chemuson.clean2d import (
    Clean2DMode,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    run_clean2d_engine,
)
from chemuson.core.model import MolGraph
from chemuson.molecular_assistant import (
    MolecularAssistant,
    MolecularAssistantRequest,
    MolecularAssistantResult,
    MolecularAssistantStatus,
    OpenAICompatibleConfig,
    OpenAICompatibleProvider,
)
from chemuson.molecular_assistant.limits import DEFAULT_PROVIDER_TIMEOUT_S


MetricRunner = Callable[..., Any]
_METRIC_FIELDS = (
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
)


def evaluate_result(
    result: MolecularAssistantResult,
    *,
    clean2d_runner: MetricRunner | None = None,
    mode: Clean2DMode | str = Clean2DMode.QUICK,
    target_bond_length: float = 42.0,
) -> dict[str, Any]:
    """Build a JSON-safe report, passing only validated successful graphs to M02."""
    if not isinstance(result, MolecularAssistantResult):
        raise TypeError("result must be a MolecularAssistantResult")
    clean2d_mode = Clean2DMode(mode)
    target = _valid_target_bond_length(target_bond_length)
    assistant_record = {
        "status": result.status.value,
        "provider_id": result.provider_id,
        "model_id": result.model_id,
        "validation_passed": result.validation_passed,
        "reason_code": result.reason_code,
        "proposed_smiles": result.proposed_smiles,
    }
    report: dict[str, Any] = {
        "schema_version": 1,
        "molecular_assistant": assistant_record,
        "clean2d": None,
    }
    if (
        result.status != MolecularAssistantStatus.SUCCESS
        or result.validation_passed is not True
        or not isinstance(result.graph, MolGraph)
        or not result.graph.atoms
    ):
        return report

    graph = result.graph
    initial_coords = _graph_coords(graph)
    before_metrics = _metrics(graph, initial_coords, target)
    # Clean2D receives a disposable clone so even a future mutating backend
    # cannot modify the validated M23 result retained by the caller.
    working_graph = copy.deepcopy(graph)
    working_snapshot = capture_clean2d_snapshot(working_graph)
    working_coords = _graph_coords(working_graph)
    runner = clean2d_runner or run_clean2d_engine
    try:
        clean_result = runner(
            working_graph,
            mode=clean2d_mode,
            target_bond_length=target,
        )
    except Exception:
        report["clean2d"] = _failed_clean2d_record(
            clean2d_mode,
            target,
            before_metrics,
            "backend-failure",
        )
        return report

    if (
        capture_clean2d_snapshot(working_graph) != working_snapshot
        or _graph_coords(working_graph) != working_coords
    ):
        report["clean2d"] = _failed_clean2d_record(
            clean2d_mode,
            target,
            before_metrics,
            "invariant-violation",
        )
        return report

    selected = getattr(clean_result, "selected", None)
    is_accepted = selected is not None and not bool(getattr(selected, "rejected", False))
    after_metrics = (
        _metrics(working_graph, getattr(selected, "coords", {}), target)
        if is_accepted
        else None
    )
    state = str(getattr(clean_result, "result_state", "failed-controlled"))
    stable_reason = str(getattr(clean_result, "stable_reason", "") or "") or None
    report["clean2d"] = {
        "mode": clean2d_mode.value,
        "target_bond_length": target,
        "state": state,
        "reason": stable_reason,
        "candidate_sources": list(getattr(clean_result, "candidate_sources", ()) or ()),
        "metrics": {"before": before_metrics, "after": after_metrics},
    }
    return report


def encode_report(report: dict[str, Any]) -> str:
    """Serialize a report with a strict JSON finite-number contract."""
    return json.dumps(report, ensure_ascii=False, sort_keys=True, allow_nan=False, indent=2)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", required=True, help="Explicit HTTP(S) OpenAI-compatible endpoint")
    parser.add_argument("--model", required=True, help="Explicit model identifier")
    parser.add_argument("--description", required=True, help="One molecular structure request")
    parser.add_argument(
        "--api-key-env",
        default="OPENAI_API_KEY",
        help="Environment variable containing an optional API key (default: OPENAI_API_KEY)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_PROVIDER_TIMEOUT_S,
        help="Finite provider timeout in seconds",
    )
    parser.add_argument(
        "--mode",
        choices=[mode.value for mode in Clean2DMode],
        default=Clean2DMode.QUICK.value,
        help="Existing Clean2D evaluation mode",
    )
    parser.add_argument(
        "--target-bond-length",
        type=float,
        default=42.0,
        help="Clean2D target bond length in canvas pixels",
    )
    args = parser.parse_args(argv)
    try:
        api_key = os.environ.get(args.api_key_env) if args.api_key_env else None
        config = OpenAICompatibleConfig(
            base_url=args.base_url,
            model=args.model,
            api_key=api_key,
            timeout_s=args.timeout,
        )
        target = _valid_target_bond_length(args.target_bond_length)
    except (TypeError, ValueError):
        parser.error("invalid provider or Clean2D configuration")

    provider = OpenAICompatibleProvider(config)
    result = MolecularAssistant(provider).generate(MolecularAssistantRequest(args.description))
    report = evaluate_result(
        result,
        mode=args.mode,
        target_bond_length=target,
    )
    print(encode_report(report))
    return 0 if result.status == MolecularAssistantStatus.SUCCESS else 1


def _metrics(graph: MolGraph, coords: dict[int, tuple[float, float]], target: float) -> dict[str, Any]:
    quality = classify_clean2d_layout_quality(
        graph,
        coords=coords,
        target_bond_length=target,
    )
    values = {name: getattr(quality, name) for name in _METRIC_FIELDS}
    return {name: _json_finite(values[name]) for name in _METRIC_FIELDS}


def _failed_clean2d_record(
    mode: Clean2DMode,
    target: float,
    before_metrics: dict[str, Any],
    reason: str,
) -> dict[str, Any]:
    return {
        "mode": mode.value,
        "target_bond_length": target,
        "state": "failed-controlled",
        "reason": reason,
        "candidate_sources": [],
        "metrics": {"before": before_metrics, "after": None},
    }


def _graph_coords(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {
        atom_id: (float(atom.x), float(atom.y))
        for atom_id, atom in graph.atoms.items()
    }


def _json_finite(value: Any) -> Any:
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    raise TypeError(f"unsupported Clean2D metric type: {type(value).__name__}")


def _valid_target_bond_length(value: float) -> float:
    if isinstance(value, bool):
        raise ValueError("target_bond_length must be finite and positive")
    try:
        target = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("target_bond_length must be finite and positive") from exc
    if not math.isfinite(target) or target <= 0:
        raise ValueError("target_bond_length must be finite and positive")
    return target


if __name__ == "__main__":
    sys.exit(main())
