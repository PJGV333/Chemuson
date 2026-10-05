from __future__ import annotations

import json
from types import SimpleNamespace

from chemuson.clean2d import Clean2DMode, run_clean2d_engine
from chemuson.core.model import MolGraph
from chemuson.molecular_assistant import (
    MolecularAssistantResult,
    MolecularAssistantStatus,
)
from tools import ai_clean2d_evaluation as evaluation


def _graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    graph.add_atom("C", 42.0, 0.0, atom_id=2)
    graph.add_bond(1, 2, bond_id=1)
    return graph


def _success(graph: MolGraph | None = None) -> MolecularAssistantResult:
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="fake-provider",
        model_id="fake-model",
        proposed_smiles="CC",
        graph=graph or _graph(),
        validation_passed=True,
    )


def _candidate_result(*, state: str = "no-op", coords=None, reason: str = ""):
    candidate = None
    if coords is not None:
        candidate = SimpleNamespace(coords=coords, rejected=False)
    return SimpleNamespace(
        selected=candidate,
        result_state=state,
        stable_reason=reason,
        candidate_sources=("current",) if candidate is not None else (),
    )


def test_success_evaluates_a_clone_and_reports_finite_before_after_metrics() -> None:
    result = _success()
    original_graph = result.graph
    original_coords = {atom_id: (atom.x, atom.y) for atom_id, atom in original_graph.atoms.items()}
    calls = []

    def fake_clean2d(graph, *, mode, target_bond_length):
        calls.append((graph, mode, target_bond_length))
        return _candidate_result(
            coords={1: (0.0, 0.0), 2: (42.0, 0.0)},
        )

    report = evaluation.evaluate_result(result, clean2d_runner=fake_clean2d)
    clean2d = report["clean2d"]
    assert len(calls) == 1
    assert calls[0][0] is not original_graph
    assert calls[0][1:] == (Clean2DMode.QUICK, 42.0)
    assert {atom_id: (atom.x, atom.y) for atom_id, atom in original_graph.atoms.items()} == original_coords
    assert clean2d["state"] == "no-op"
    assert clean2d["reason"] is None
    assert set(clean2d["metrics"]["before"]) == set(evaluation._METRIC_FIELDS)
    assert set(clean2d["metrics"]["after"]) == set(evaluation._METRIC_FIELDS)
    assert clean2d["metrics"]["before"]["min_nonbonded_distance"] is None
    assert clean2d["metrics"]["before"]["min_ring_degeneracy"] is None
    decoded = json.loads(evaluation.encode_report(report))
    assert decoded["molecular_assistant"]["provider_id"] == "fake-provider"
    assert "description" not in decoded["molecular_assistant"]


def test_failed_molecular_result_never_reaches_clean2d() -> None:
    result = MolecularAssistantResult(
        status=MolecularAssistantStatus.INVALID_STRUCTURE,
        provider_id="fake-provider",
        model_id="fake-model",
        proposed_smiles="invalid",
        graph=None,
        validation_passed=False,
        reason_code="invalid_smiles",
    )

    def unexpected_runner(*_args, **_kwargs):
        raise AssertionError("Clean2D must not receive a failed M23 result")

    report = evaluation.evaluate_result(result, clean2d_runner=unexpected_runner)
    assert report["molecular_assistant"]["status"] == "invalid_structure"
    assert report["molecular_assistant"]["reason_code"] == "invalid_smiles"
    assert report["clean2d"] is None
    assert "metrics" not in evaluation.encode_report(report)


def test_no_accepted_candidate_keeps_engine_state_and_omits_after_metrics() -> None:
    report = evaluation.evaluate_result(
        _success(),
        clean2d_runner=lambda *_args, **_kwargs: _candidate_result(
            state="failed-controlled",
            reason="backend-failure",
        ),
    )

    assert report["clean2d"]["state"] == "failed-controlled"
    assert report["clean2d"]["reason"] == "backend-failure"
    assert report["clean2d"]["metrics"]["before"]
    assert report["clean2d"]["metrics"]["after"] is None


def test_metrics_do_not_override_clean2d_state() -> None:
    report = evaluation.evaluate_result(
        _success(),
        clean2d_runner=lambda *_args, **_kwargs: _candidate_result(
            state="no-op",
            coords={1: (0.0, 0.0), 2: (500.0, 0.0)},
        ),
    )

    assert report["clean2d"]["state"] == "no-op"
    assert report["clean2d"]["metrics"]["after"]["quality_class"] == "needs_rebuild"


def test_mutating_engine_is_isolated_and_reported_as_invariant_violation() -> None:
    result = _success()
    original_coords = {atom_id: (atom.x, atom.y) for atom_id, atom in result.graph.atoms.items()}

    def mutating_runner(graph, **_kwargs):
        graph.atoms[1].x = 999.0
        return _candidate_result(coords={1: (0.0, 0.0), 2: (42.0, 0.0)})

    report = evaluation.evaluate_result(result, clean2d_runner=mutating_runner)
    assert report["clean2d"]["state"] == "failed-controlled"
    assert report["clean2d"]["reason"] == "invariant-violation"
    assert {atom_id: (atom.x, atom.y) for atom_id, atom in result.graph.atoms.items()} == original_coords


def test_real_clean2d_engine_evaluates_a_small_graph_with_bounded_backend() -> None:
    def bounded_runner(graph, **kwargs):
        return run_clean2d_engine(graph, rdkit_timeout_s=0.25, seed=1, **kwargs)

    report = evaluation.evaluate_result(_success(), clean2d_runner=bounded_runner)
    clean2d = report["clean2d"]
    assert clean2d["state"] in {"applied", "no-op", "preserve-only", "failed-controlled"}
    assert clean2d["metrics"]["before"]["quality_class"]
    if clean2d["state"] != "failed-controlled":
        assert clean2d["metrics"]["after"] is not None
    json.dumps(report, allow_nan=False)


def test_cli_requires_explicit_provider_and_request_fields() -> None:
    try:
        evaluation.main([])
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("CLI should reject missing endpoint, model, and description")
