from __future__ import annotations

import json

from chemuson.clean2d import (
    describe_clean2d_rigid_systems,
    generate_clean2d_candidates,
    plan_clean2d_block_assembly,
    run_clean2d_engine,
)
from tests.clean2d_regression.campaign4_cases import CAMPAIGN4_CASES, CASES_BY_NAME
from tests.clean2d_regression.cases import get_regression_cases


_EXPECTED_FAMILIES = {
    "regular_benzene_control": "monocycle",
    "simple_monocycle_control": "monocycle",
    "fused_aromatic_bicyclic": "fused",
    "fused_non_aromatic_bicyclic": "fused",
    "spiro_bicyclic": "spiro",
    "bridged_bicyclic": "bridged",
    "larger_polycyclic": "fused",
    "fused_one_substituent": "fused",
    "fused_multiple_substituents": "fused",
    "spiro_substituent": "spiro",
    "two_rigid_blocks_linker": "monocycle",
    "congested_ring_attachments": "monocycle",
}


def test_campaign4_rigid_descriptor_is_complete_deterministic_and_json_safe() -> None:
    for case in CAMPAIGN4_CASES:
        graph = case.builder()
        first = describe_clean2d_rigid_systems(graph)
        second = describe_clean2d_rigid_systems(graph)
        assert first == second, case.name
        assert first["systems"], case.name
        assert any(system["family"] == _EXPECTED_FAMILIES[case.name] for system in first["systems"]), case.name
        for system in first["systems"]:
            assert set(system) >= {
                "id",
                "family",
                "atom_ids",
                "bond_ids",
                "ring_ids",
                "shared_atom_ids",
                "shared_bond_ids",
                "attachment_atom_ids",
                "external_neighbor_ids",
                "centroid",
                "principal_orientation",
                "attachment_vectors",
                "external_substituent_count",
                "local_congestion",
            }
        json.dumps(first, allow_nan=False, sort_keys=True)


def test_campaign4_emits_internal_rigid_multiring_candidate_with_metadata() -> None:
    for case_name in ("fused_aromatic_bicyclic", "spiro_bicyclic", "bridged_bicyclic", "congested_ring_attachments"):
        graph = CASES_BY_NAME[case_name].builder()
        candidates = generate_clean2d_candidates(graph, mode="quick", target_bond_length=40.0)
        candidate = next(item for item in candidates if item.source == "rigid_multiring_layout")
        assert candidate.metadata["strategy"] == "rigid_multiring_layout"
        assert candidate.metadata["rigid_system_count"] >= 1
        assert candidate.metadata["rigid_system_types"]
        assert candidate.metadata["affected_atom_ids"]
        assert candidate.metadata["attachment_atoms"]
        assert set(candidate.metadata["hard_gate_checks"]) >= {
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
        json.dumps(candidate.metadata, allow_nan=False, sort_keys=True)


def test_campaign4_controls_do_not_force_internal_candidate() -> None:
    for case_name in ("regular_benzene_control", "simple_monocycle_control"):
        candidates = generate_clean2d_candidates(CASES_BY_NAME[case_name].builder(), mode="quick", target_bond_length=40.0)
        assert not any(candidate.source == "rigid_multiring_layout" for candidate in candidates)


def test_campaign4_multiple_rigid_blocks_emit_local_candidate_without_replanning() -> None:
    candidates = generate_clean2d_candidates(
        CASES_BY_NAME["two_rigid_blocks_linker"].builder(),
        mode="quick",
        target_bond_length=40.0,
    )
    candidate = next(item for item in candidates if item.source == "rigid_multiring_layout")
    assert candidate.metadata["assembly_plan"]["requires_global_assembly"] is True
    assert candidate.metadata["layout_mode"] == "attachment_orientation"


def test_campaign4_candidate_ordering_and_result_are_deterministic() -> None:
    graph = CASES_BY_NAME["fused_multiple_substituents"].builder()
    first = generate_clean2d_candidates(graph, mode="quick", target_bond_length=40.0)
    second = generate_clean2d_candidates(graph, mode="quick", target_bond_length=40.0)
    assert [(item.source, item.rejected, item.rejection_reason) for item in first] == [
        (item.source, item.rejected, item.rejection_reason) for item in second
    ]
    first_result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    second_result = run_clean2d_engine(CASES_BY_NAME["fused_multiple_substituents"].builder(), mode="quick", target_bond_length=40.0)
    assert first_result.result_state == second_result.result_state
    assert first_result.selected is not None and second_result.selected is not None
    assert first_result.selected.source == second_result.selected.source


def test_campaign4_preserves_campaign3_global_placement_controls() -> None:
    cases = {case.name: case for case in get_regression_cases()}
    for name in ("multiblock_triphenyl_like", "multiblock_branched"):
        result = run_clean2d_engine(cases[name].builder(), mode="quick", target_bond_length=40.0)
        assert result.selected is not None
        assert result.selected.source == "global_block_placement"
        assert result.result_state != "preserve-only"


def test_campaign4_consumes_campaign3_plan_without_reassigning_blocks() -> None:
    plan = plan_clean2d_block_assembly(CASES_BY_NAME["two_rigid_blocks_linker"].builder())
    assert plan["requires_global_assembly"] is True
    assert plan["anchor_block_id"] is not None
