from __future__ import annotations

import pytest

from chemuson.clean2d import (
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    run_clean2d_engine,
)
from tests.clean2d_regression.assertions import execute_case
from tests.clean2d_regression.cases import get_regression_cases


SIMPLE_AROMATIC_CASES = (
    "aromatic_benzene_regular",
    "heteroaromatic_pyridine_like",
    "heteroaromatic_furan_like",
    "heteroaromatic_thiophene_like",
)


def _case(name: str):
    return next(case for case in get_regression_cases() if case.name == name)


def _coords(graph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


@pytest.mark.parametrize("case_name", SIMPLE_AROMATIC_CASES)
def test_simple_isolated_aromatic_ring_improves_measurably(case_name: str) -> None:
    case = _case(case_name)
    graph = case.builder()
    before = _coords(graph)
    before_quality = classify_clean2d_layout_quality(graph, coords=before, target_bond_length=40.0)

    result = run_clean2d_engine(graph, mode=case.mode, target_bond_length=40.0)

    assert result.result_state == "applied"
    assert result.selected is not None
    assert result.selected.source == "simple_aromatic_template"
    after_quality = classify_clean2d_layout_quality(graph, coords=result.selected.coords, target_bond_length=40.0)
    assert after_quality.length_rms_error < before_quality.length_rms_error
    assert after_quality.length_max_error < before_quality.length_max_error
    assert after_quality.visual_score < before_quality.visual_score
    assert after_quality.crossings <= before_quality.crossings
    assert after_quality.min_ring_degeneracy >= before_quality.min_ring_degeneracy - 1e-9


def test_simple_aromatic_improvement_preserves_chemical_identity() -> None:
    case = _case("aromatic_benzene_regular")
    graph = case.builder()
    before = _coords(graph)
    snapshot = capture_clean2d_snapshot(graph)

    result = run_clean2d_engine(graph, mode=case.mode, target_bond_length=40.0)

    assert result.selected is not None
    assert_clean2d_invariants(snapshot, graph, before, result.selected.coords)


def test_simple_aromatic_candidate_does_not_route_complex_policy_guards() -> None:
    guarded_cases = [case for case in get_regression_cases() if case.has_tag("complex_policy_guard")]
    assert guarded_cases

    for case in guarded_cases:
        snapshot = execute_case(case)
        assert "simple_aromatic_template" not in snapshot["result"]["candidate_sources"]
        assert snapshot["result"]["state"] in case.expected_states
