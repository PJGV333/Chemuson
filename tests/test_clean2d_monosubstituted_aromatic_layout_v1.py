from __future__ import annotations

import math

import pytest

from chemuson.clean2d import (
    Clean2DMode,
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
    classify_clean2d_layout_quality,
    generate_clean2d_candidates,
    run_clean2d_engine,
)
from tests.clean2d_regression.assertions import execute_case
from tests.clean2d_regression.cases import get_regression_cases


MONOSUBSTITUTED_CASES = (
    "aromatic_toluene_like",
    "aromatic_phenol_like",
    "aromatic_aniline_like",
)


def _case(name: str):
    return next(case for case in get_regression_cases() if case.name == name)


def _coords(graph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _ring_and_substituent(graph) -> tuple[set[int], int, int]:
    aromatic_atoms = {
        atom_id
        for bond in graph.bonds.values()
        if bool(getattr(bond, "is_aromatic", False))
        for atom_id in (bond.a1_id, bond.a2_id)
    }
    substituents = sorted(set(graph.atoms) - aromatic_atoms)
    assert len(substituents) == 1
    substituent_id = substituents[0]
    anchor_id = next(
        bond.a1_id if bond.a2_id == substituent_id else bond.a2_id
        for bond in graph.bonds.values()
        if substituent_id in {bond.a1_id, bond.a2_id}
    )
    return aromatic_atoms, anchor_id, substituent_id


def _dot(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1]


@pytest.mark.parametrize("case_name", MONOSUBSTITUTED_CASES)
def test_monosubstituted_aromatic_improves_measurably(case_name: str) -> None:
    case = _case(case_name)
    graph = case.builder()
    before = _coords(graph)
    before_quality = classify_clean2d_layout_quality(graph, coords=before, target_bond_length=40.0)

    result = run_clean2d_engine(graph, mode=case.mode, target_bond_length=40.0)

    assert result.result_state == "applied"
    assert result.selected is not None
    assert result.selected.source == "monosubstituted_aromatic_template"
    after_quality = classify_clean2d_layout_quality(graph, coords=result.selected.coords, target_bond_length=40.0)
    assert after_quality.visual_score < before_quality.visual_score
    assert after_quality.length_rms_error < before_quality.length_rms_error
    assert after_quality.length_max_error < before_quality.length_max_error
    assert after_quality.angle_rms_deviation <= before_quality.angle_rms_deviation
    assert after_quality.crossings <= before_quality.crossings
    assert after_quality.min_ring_degeneracy >= before_quality.min_ring_degeneracy - 1e-9
    assert after_quality.min_nonbonded_distance > 40.0 * 0.50


def test_monosubstituted_aromatic_substituent_points_outward() -> None:
    case = _case("aromatic_toluene_like")
    graph = case.builder()
    ring, anchor_id, substituent_id = _ring_and_substituent(graph)

    result = run_clean2d_engine(graph, mode=case.mode, target_bond_length=40.0)

    assert result.selected is not None
    coords = result.selected.coords
    center = (
        sum(coords[atom_id][0] for atom_id in ring) / len(ring),
        sum(coords[atom_id][1] for atom_id in ring) / len(ring),
    )
    outward = (coords[anchor_id][0] - center[0], coords[anchor_id][1] - center[1])
    substituent = (
        coords[substituent_id][0] - coords[anchor_id][0],
        coords[substituent_id][1] - coords[anchor_id][1],
    )

    assert _dot(outward, substituent) > 0.0
    assert math.hypot(*substituent) == pytest.approx(40.0, rel=0.08)


def test_monosubstituted_aromatic_preserves_chemical_identity() -> None:
    case = _case("aromatic_phenol_like")
    graph = case.builder()
    before = _coords(graph)
    snapshot = capture_clean2d_snapshot(graph)

    result = run_clean2d_engine(graph, mode=case.mode, target_bond_length=40.0)

    assert result.selected is not None
    assert_clean2d_invariants(snapshot, graph, before, result.selected.coords)


def test_monosubstituted_aromatic_template_is_not_a_propose_candidate() -> None:
    case = _case("aromatic_toluene_like")
    graph = case.builder()

    candidates = generate_clean2d_candidates(graph, mode=Clean2DMode.PROPOSE, target_bond_length=40.0)

    assert "monosubstituted_aromatic_template" not in {candidate.source for candidate in candidates}


def test_monosubstituted_aromatic_template_excludes_disubstituted_and_fused_cases() -> None:
    excluded = (
        "aromatic_para_disubstituted_like",
        "fused_aromatic_current_baseline",
        "fused_aromatic_anthracene_like",
    )
    for case_name in excluded:
        case = _case(case_name)
        graph = case.builder()
        candidates = generate_clean2d_candidates(graph, mode=case.mode, target_bond_length=40.0)
        assert "monosubstituted_aromatic_template" not in {candidate.source for candidate in candidates}


def test_monosubstituted_aromatic_template_does_not_route_complex_policy_guards() -> None:
    guarded_cases = [case for case in get_regression_cases() if case.has_tag("complex_policy_guard")]
    assert guarded_cases

    for case in guarded_cases:
        snapshot = execute_case(case)
        assert "monosubstituted_aromatic_template" not in snapshot["result"]["candidate_sources"]


def test_simple_aromatic_template_still_handles_isolated_aromatics() -> None:
    for case_name in (
        "aromatic_benzene_regular",
        "heteroaromatic_pyridine_like",
        "heteroaromatic_furan_like",
        "heteroaromatic_thiophene_like",
    ):
        case = _case(case_name)
        result = run_clean2d_engine(case.builder(), mode=case.mode, target_bond_length=40.0)
        assert result.selected is not None
        assert result.selected.source == "simple_aromatic_template"
