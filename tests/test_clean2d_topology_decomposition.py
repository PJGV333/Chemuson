from __future__ import annotations

import json
import math

from chemuson.clean2d.complex_policy import describe_clean2d_topology
from chemuson.core.layers import build_multilayer_chemical_graph
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
    assert json.dumps(first, sort_keys=True, separators=(",", ":")) == json.dumps(
        second,
        sort_keys=True,
        separators=(",", ":"),
    )
    assert first["block_kind_counts"] == dict(sorted(first["block_kind_counts"].items()))
    assert first["components"] == sorted(first["components"], key=lambda component: tuple(component))


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


def _add_hexagon(graph: MolGraph, cx: float, cy: float, radius: float = 24.0) -> list[int]:
    atoms = []
    for index in range(6):
        angle = math.radians(60.0 * index)
        atoms.append(graph.add_atom("C", cx + math.cos(angle) * radius, cy + math.sin(angle) * radius).id)
    for index, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(index + 1) % len(atoms)], is_aromatic=True)
    return atoms
