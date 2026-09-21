from __future__ import annotations

import json
from typing import cast

import pytest

from chemuson.clean2d import generate_clean2d_candidates, plan_clean2d_block_assembly
from tests.clean2d_regression.assertions import execute_case
from tests.clean2d_regression.cases import get_regression_cases


_CASES = {case.name: case for case in get_regression_cases()}


@pytest.mark.parametrize(
    ("case_name", "requires_global_assembly"),
    (
        ("multiblock_biphenyl_like", True),
        ("multiblock_diphenyl_ether_like", True),
        ("multiblock_triphenyl_like", True),
        ("aromatic_benzene_regular", False),
        ("acyclic_butane_stretched", False),
    ),
)
def test_medium_assembly_plan_is_topology_derived_and_json_safe(
    case_name: str,
    requires_global_assembly: bool,
) -> None:
    graph = _CASES[case_name].builder()
    first = plan_clean2d_block_assembly(graph)
    second = plan_clean2d_block_assembly(graph)

    assert first == second
    assert first["requires_global_assembly"] is requires_global_assembly
    ordered_block_ids = cast(list[int], first["ordered_block_ids"])
    assert len(ordered_block_ids) == len(set(ordered_block_ids))
    ordered_connector_ids = cast(list[int], first["ordered_connector_ids"])
    assert len(ordered_connector_ids) == len(set(ordered_connector_ids))
    assert isinstance(first["anchor_block_id"], int) if requires_global_assembly else first["anchor_block_id"] is None
    json.dumps(first, allow_nan=False, sort_keys=True)


def test_medium_block_candidate_records_global_assembly_plan() -> None:
    graph = _CASES["multiblock_diphenyl_ether_like"].builder()
    candidates = generate_clean2d_candidates(graph, mode="publication", target_bond_length=40.0)
    block_candidate = next(candidate for candidate in candidates if candidate.source in {"block_constraints", "block_layout"})

    plan = block_candidate.metadata["global_assembly_plan"]
    assert plan["requires_global_assembly"] is True
    assert plan["ordered_connector_ids"]
    json.dumps(plan, allow_nan=False, sort_keys=True)


@pytest.mark.parametrize(
    "case_name",
    ("multiblock_biphenyl_like", "multiblock_diphenyl_ether_like", "aromatic_benzene_regular", "acyclic_butane_stretched"),
)
def test_medium_assembly_regression_controls_preserve_contract(case_name: str) -> None:
    record = execute_case(_CASES[case_name])
    assert record["identity"]["atom_ids"]
    assert record["identity"]["bond_ids"]
    assert record["result"]["state"] in _CASES[case_name].expected_states
    if record["metrics"]["after"] is not None and "multiblock" in case_name:
        assert record["metrics"]["after"]["crossings"] == 0
