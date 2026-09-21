from __future__ import annotations

import json
import math
from itertools import pairwise

from chemuson.clean2d.complex_policy import (
    classify_clean2d_complexity,
    describe_clean2d_topology,
)
from chemuson.clean2d.engine import run_clean2d_engine
from chemuson.core.layers import BlockKind, build_multilayer_chemical_graph
from chemuson.core.model import MolGraph


def test_topology_summary_exposes_existing_blocks_and_connectors() -> None:
    graph = _two_ring_linker_graph()
    layers = build_multilayer_chemical_graph(graph)

    summary = describe_clean2d_topology(layers)

    assert summary["atom_count"] == len(graph.atoms)
    assert summary["bond_count"] == len(graph.bonds)
    assert summary["component_count"] == 1
    assert summary["ring_count"] >= 2
    assert summary["block_count"] == len(layers.block_graph.blocks)
    assert summary["connector_count"] == len(layers.block_graph.edges)
    assert summary["blocks"]
    assert all({"id", "kind", "atom_ids", "anchor_atom_ids", "motif_ids"} <= block.keys() for block in summary["blocks"])
    assert all({"id", "kind", "block_ids", "atom_ids", "weight"} <= edge.keys() for edge in summary["connectors"])


def test_topology_summary_is_deterministic_and_json_safe() -> None:
    graph = _two_ring_linker_graph()
    layers = build_multilayer_chemical_graph(graph)

    first = describe_clean2d_topology(layers)
    second = describe_clean2d_topology(layers)

    assert first == second
    assert json.dumps(first, allow_nan=False, sort_keys=True, separators=(",", ":")) == json.dumps(
        second,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    assert first["block_kind_counts"] == dict(sorted(first["block_kind_counts"].items()))
    assert first["components"] == sorted(first["components"], key=lambda component: tuple(component))


def test_topology_matrix_acyclic_chain() -> None:
    summary = _summary(_chain_graph(5))
    assert summary["ring_count"] == 0
    assert summary["component_count"] == 1
    assert summary["articulation_atom_ids"] == [2, 3, 4]
    assert summary["branch_point_atom_ids"] == []


def test_topology_matrix_branched_tree() -> None:
    graph = _chain_graph(4)
    branch = graph.add_atom("C", 0.0, 20.0).id
    graph.add_bond(2, branch, order=1)
    summary = _summary(graph)
    assert summary["ring_count"] == 0
    assert summary["branch_point_atom_ids"] == [2]
    assert 2 in summary["articulation_atom_ids"]


def test_topology_matrix_monocycle() -> None:
    summary = _summary(_cycle_graph(6))
    assert summary["ring_count"] == 1
    assert summary["ring_systems"][0]["kind"] == "monocycle"
    assert summary["fused_systems"] == []
    assert summary["spiro_systems"] == []


def test_topology_matrix_fused_bicyclic() -> None:
    summary = _summary(_fused_graph())
    assert summary["ring_count"] >= 2
    assert summary["fused_systems"]
    assert any(block["kind"] == BlockKind.FUSED_SYSTEM.value for block in summary["blocks"])


def test_topology_matrix_spiro_system() -> None:
    summary = _summary(_spiro_graph())
    assert summary["spiro_systems"]
    assert summary["spiro_systems"][0]["shared_atom_ids"]


def test_topology_matrix_bridged_system() -> None:
    summary = _summary(_bridged_graph())
    assert summary["bridged_systems"]
    assert summary["bridged_systems"][0]["bridgehead_atom_ids"]


def test_topology_matrix_two_rings_and_flexible_linker() -> None:
    summary = _summary(_two_ring_linker_graph())
    assert summary["flexible_connectors"]
    assert summary["attachment_atom_ids"]
    assert summary["block_adjacency"]


def test_topology_matrix_multiblock_structure() -> None:
    graph = _two_ring_linker_graph()
    third = _add_hexagon(graph, 360.0, 0.0)
    graph.add_bond(15, third[0], order=1)
    summary = _summary(graph)
    assert summary["block_count"] >= 3
    assert len(summary["block_adjacency"]) >= 2
    assert len(summary["flexible_connectors"]) >= 2


def test_topology_matrix_macrocycle() -> None:
    summary = _summary(_cycle_graph(12))
    assert summary["macrocycle_blocks"]
    assert summary["ring_count"] >= 1


def test_topology_summary_handles_disconnected_components_and_selection() -> None:
    graph = _chain_graph(3)
    second = graph.add_atom("C", 100.0, 0.0).id
    third = graph.add_atom("C", 140.0, 0.0).id
    graph.add_bond(second, third, order=1)
    summary = _summary(graph)
    assert summary["component_count"] == 2
    selected = build_multilayer_chemical_graph(graph, {1, 2, 3})
    selected_summary = describe_clean2d_topology(selected)
    assert selected_summary["atom_ids"] == [1, 2, 3]
    assert selected_summary["component_count"] == 1


def test_topology_summary_is_strict_json_safe_and_observational() -> None:
    graph = _two_ring_linker_graph()
    layers = build_multilayer_chemical_graph(graph)
    before_atoms = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    before_bonds = {bond_id: (bond.a1_id, bond.a2_id, bond.order) for bond_id, bond in graph.bonds.items()}
    layers.block_graph.blocks[0].metadata["nonfinite"] = float("nan")

    summary = describe_clean2d_topology(layers)

    json.dumps(summary, allow_nan=False, sort_keys=True)
    assert summary["blocks"][0]["metadata"]["nonfinite"] is None
    assert before_atoms == {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    assert before_bonds == {bond_id: (bond.a1_id, bond.a2_id, bond.order) for bond_id, bond in graph.bonds.items()}


def test_topology_evidence_does_not_change_policy_or_candidate_outcome() -> None:
    graph = _two_ring_linker_graph()
    before_profile = classify_clean2d_complexity(graph)
    before_result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    describe_clean2d_topology(build_multilayer_chemical_graph(graph))

    after_profile = classify_clean2d_complexity(graph)
    after_result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)
    assert before_profile.policy_evidence == after_profile.policy_evidence
    assert _result_signature(before_result) == _result_signature(after_result)


def _result_signature(result: object) -> tuple[object, ...]:
    selected = getattr(result, "selected", None)
    candidates = tuple(getattr(result, "candidates", ()) or ())
    rejected = tuple(getattr(result, "rejected", ()) or ())
    return (
        bool(getattr(result, "ok", False)),
        getattr(selected, "source", None),
        tuple(candidate.source for candidate in candidates),
        tuple(candidate.source for candidate in rejected),
    )


def _two_ring_linker_graph() -> MolGraph:
    graph = MolGraph()
    left = _add_hexagon(graph, 0.0, 0.0)
    right = _add_hexagon(graph, 180.0, 0.0)
    linker_a = graph.add_atom("C", 54.0, 0.0).id
    linker_b = graph.add_atom("C", 126.0, 0.0).id
    graph.add_bond(left[0], linker_a, order=1)
    graph.add_bond(linker_a, linker_b, order=1)
    graph.add_bond(linker_b, right[3], order=1)
    return graph


def _summary(graph: MolGraph) -> dict[str, object]:
    return describe_clean2d_topology(build_multilayer_chemical_graph(graph))


def _chain_graph(atom_count: int) -> MolGraph:
    graph = MolGraph()
    atoms = [graph.add_atom("C", float(index) * 40.0, 0.0).id for index in range(1, atom_count + 1)]
    for left, right in pairwise(atoms):
        graph.add_bond(left, right, order=1)
    return graph


def _cycle_graph(atom_count: int) -> MolGraph:
    graph = MolGraph()
    atoms = _add_cycle_atoms(graph, atom_count, 0.0, 0.0)
    _close_cycle(graph, atoms)
    return graph


def _fused_graph() -> MolGraph:
    graph = MolGraph()
    left = _add_cycle_atoms(graph, 6, 0.0, 0.0)
    right = [left[2], left[3], *_add_cycle_atoms(graph, 4, 48.0, 0.0)]
    _close_cycle(graph, left)
    _close_cycle(graph, right)
    return graph


def _spiro_graph() -> MolGraph:
    graph = MolGraph()
    center = graph.add_atom("C", 0.0, 0.0).id
    first = [center, *_add_cycle_atoms(graph, 4, 20.0, 0.0)]
    second = [center, *_add_cycle_atoms(graph, 4, -20.0, 0.0)]
    _close_cycle(graph, first)
    _close_cycle(graph, second)
    return graph


def _bridged_graph() -> MolGraph:
    graph = MolGraph()
    bridgeheads = [graph.add_atom("C", 0.0, 0.0).id, graph.add_atom("C", 120.0, 0.0).id]
    paths = []
    for y, length in ((40.0, 2), (-40.0, 2), (0.0, 1)):
        middle = [graph.add_atom("C", 40.0 + index * 30.0, y).id for index in range(length)]
        path = [bridgeheads[0], *middle, bridgeheads[1]]
        paths.append(path)
        for left, right in pairwise(path):
            graph.add_bond(left, right, order=1)
    return graph


def _add_cycle_atoms(graph: MolGraph, atom_count: int, cx: float, cy: float) -> list[int]:
    return [
        graph.add_atom(
            "C",
            cx + math.cos(math.radians(360.0 * index / atom_count)) * 24.0,
            cy + math.sin(math.radians(360.0 * index / atom_count)) * 24.0,
        ).id
        for index in range(atom_count)
    ]


def _close_cycle(graph: MolGraph, atoms: list[int]) -> None:
    for left, right in zip(atoms, atoms[1:] + atoms[:1]):
        graph.add_bond(left, right, order=1)


def _add_hexagon(graph: MolGraph, cx: float, cy: float, radius: float = 24.0) -> list[int]:
    atoms = []
    for index in range(6):
        angle = math.radians(60.0 * index)
        atoms.append(graph.add_atom("C", cx + math.cos(angle) * radius, cy + math.sin(angle) * radius).id)
    for index, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(index + 1) % len(atoms)], is_aromatic=True)
    return atoms
