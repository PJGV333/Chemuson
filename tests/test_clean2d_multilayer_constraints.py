from __future__ import annotations

import math
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import generate_clean2d_candidates, run_clean2d_engine
from chemuson.core.model import BondStyle, MolGraph


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def test_clean2d_moves_motif_blocks_for_coordination_without_valence_change() -> None:
    graph = MolGraph()
    metal = graph.add_atom("Fe", 0.0, 0.0, is_coordination_center=True).id
    n = graph.add_atom("N", 200.0, 0.0).id
    c = graph.add_atom("C", 240.0, 0.0).id
    graph.add_bond(n, c, order=1)
    graph.add_bond(metal, n, style=BondStyle.COORDINATION, donor_atom_id=n, length_px=48.0)

    before_valence = graph.bond_order_sum(metal)
    result = run_clean2d_engine(graph, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert result.selected.source == "motif_constraints"
    assert graph.bond_order_sum(metal) == before_valence == 0.0
    after = result.selected.coords
    assert _distance(after[metal], after[n]) < 80.0
    assert math.isclose(_distance(after[n], after[c]), 40.0, rel_tol=0.25)


def _add_hexagon(graph: MolGraph, cx: float, cy: float, radius: float = 24.0) -> list[int]:
    atoms = []
    for idx in range(6):
        angle = math.radians(60.0 * idx)
        atoms.append(graph.add_atom("C", cx + math.cos(angle) * radius, cy + math.sin(angle) * radius).id)
    for idx, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(idx + 1) % len(atoms)], is_aromatic=True)
    return atoms


def test_clean2d_generates_internal_block_transform_across_linker() -> None:
    graph = MolGraph()
    left_ring = _add_hexagon(graph, 0.0, 0.0)
    right_ring = _add_hexagon(graph, 180.0, 0.0)
    linker_a = graph.add_atom("C", 54.0, 0.0).id
    linker_b = graph.add_atom("C", 126.0, 0.0).id
    graph.add_bond(left_ring[0], linker_a, order=1)
    graph.add_bond(linker_a, linker_b, order=1)
    graph.add_bond(linker_b, right_ring[3], order=1)
    graph.add_bond(left_ring[2], right_ring[1], style=BondStyle.INTERACTION, length_px=70.0)

    before = {
        atom_id: (atom.x, atom.y)
        for atom_id, atom in graph.atoms.items()
    }
    before_distance = _distance(before[left_ring[2]], before[right_ring[1]])
    candidates = generate_clean2d_candidates(graph, mode="publication", target_bond_length=40.0)
    block_candidate = next(candidate for candidate in candidates if candidate.source == "block_constraints")

    assert block_candidate.metadata["block_operation_count"] >= 1
    assert block_candidate.metadata["interaction_constraint_error"] < block_candidate.metadata["interaction_constraint_error_before"]
    assert _distance(block_candidate.coords[left_ring[2]], block_candidate.coords[right_ring[1]]) < before_distance
