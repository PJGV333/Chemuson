from __future__ import annotations

import json

from chemuson.clean2d import (
    describe_clean2d_rigid_systems,
    generate_clean2d_candidates,
    run_clean2d_engine,
)
from tests.clean2d_regression.campaign4_cases import (
    CASES_BY_NAME,
    _three_fused_ring_polycyclic,
)


_REQUIRED_GATES = {
    "finite_coordinates",
    "selection_invariants",
    "stereo_signature",
    "no_new_crossings",
    "collision_safety",
    "ring_degeneracy",
    "bond_length_sanity",
    "bounding_box_sanity",
    "displacement_budget",
}


def test_corrective_candidate_persists_boolean_gates_without_ranking_states() -> None:
    candidate = next(
        item
        for item in generate_clean2d_candidates(
            CASES_BY_NAME["fused_non_aromatic_bicyclic"].builder(),
            mode="quick",
            target_bond_length=40.0,
        )
        if item.source == "rigid_multiring_layout"
    )
    gates = candidate.metadata["hard_gate_checks"]
    assert isinstance(gates, dict)
    assert set(gates) >= _REQUIRED_GATES
    assert all(isinstance(value, bool) for value in gates.values())
    assert candidate.metadata["hard_gates_passed"] is all(gates.values())
    assert "accepted_by_engine" not in candidate.metadata
    assert "selected" not in candidate.metadata
    json.dumps(candidate.metadata, allow_nan=False, sort_keys=True)
    result = run_clean2d_engine(CASES_BY_NAME["fused_non_aromatic_bicyclic"].builder(), mode="quick", target_bond_length=40.0)
    evaluated = next(item for item in result.candidates if item.source == "rigid_multiring_layout")
    assert evaluated.metadata["accepted_by_engine"] is True
    assert evaluated.metadata["selected"] is False


def test_corrective_descriptor_reports_true_polycyclic_and_multiple_rigid_blocks() -> None:
    polycyclic = describe_clean2d_rigid_systems(_three_fused_ring_polycyclic())
    assert any(system["family"] == "polycyclic" for system in polycyclic["systems"])
    assert max(system["ring_count"] for system in polycyclic["systems"]) >= 3

    multiple = describe_clean2d_rigid_systems(CASES_BY_NAME["two_rigid_blocks_linker"].builder())
    assert multiple["multiple_rigid_blocks"] is True
    assert multiple["rigid_system_count"] >= 2


def test_corrective_spiro_target_is_selected_after_real_sector_improvement() -> None:
    result = run_clean2d_engine(
        CASES_BY_NAME["spiro_bicyclic"].builder(),
        mode="quick",
        target_bond_length=40.0,
    )
    assert result.selected is not None
    assert result.selected.source == "rigid_multiring_layout"
    candidate = result.selected
    assert candidate.metadata["hard_gates_passed"] is True
    assert candidate.metadata["metrics_after"]["crossings"] < candidate.metadata["metrics_before"]["crossings"]


def test_corrective_fused_substitution_rejects_unsafe_global_motion() -> None:
    for name in ("fused_one_substituent", "fused_multiple_substituents"):
        candidate = next(
            item
            for item in generate_clean2d_candidates(
                CASES_BY_NAME[name].builder(),
                mode="quick",
                target_bond_length=40.0,
            )
            if item.source == "rigid_multiring_layout"
        )
        assert candidate.rejected is True
