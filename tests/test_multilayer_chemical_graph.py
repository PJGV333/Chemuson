from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.layers import (
    InteractionKind,
    LayoutConstraintKind,
    MotifKind,
    build_multilayer_chemical_graph,
)
from chemuson.core.model import BondStyle, MolGraph


def test_multilayer_graph_separates_covalence_from_interactions() -> None:
    graph = MolGraph()
    metal = graph.add_atom("Fe", 0.0, 0.0, is_coordination_center=True).id
    n = graph.add_atom("N", 120.0, 0.0).id
    c = graph.add_atom("C", 160.0, 0.0).id
    graph.add_bond(n, c, order=1)
    graph.add_bond(metal, n, style=BondStyle.COORDINATION, donor_atom_id=n)

    before_valence = graph.bond_order_sum(metal)
    layers = build_multilayer_chemical_graph(graph)

    assert graph.bond_order_sum(metal) == before_valence == 0.0
    assert [edge.kind for edge in layers.interaction_graph.edges] == [InteractionKind.COORDINATION]
    assert any(motif.kind == MotifKind.LIGAND and n in motif.atom_ids for motif in layers.motif_graph.motifs)
    assert any(
        constraint.kind == LayoutConstraintKind.INTERACTION_DISTANCE
        and constraint.atom_ids == (metal, n)
        for constraint in layers.layout_constraint_graph.constraints
    )


def test_multilayer_graph_detects_ring_centroid_and_rigid_block() -> None:
    graph = MolGraph()
    atoms = [graph.add_atom("C", float(idx), 0.0).id for idx in range(6)]
    for idx, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(idx + 1) % len(atoms)], is_aromatic=True)

    layers = build_multilayer_chemical_graph(graph)
    kinds = [motif.kind for motif in layers.motif_graph.motifs]

    assert MotifKind.RING in kinds
    assert MotifKind.CENTROID in kinds
    assert MotifKind.RIGID_BLOCK in kinds
    assert any(
        constraint.kind == LayoutConstraintKind.MOTIF_RIGID
        for constraint in layers.layout_constraint_graph.constraints
    )


def test_explicit_non_covalent_interaction_kinds_become_layout_constraints() -> None:
    graph = MolGraph()
    kinds = [
        InteractionKind.HYDROGEN_BOND,
        InteractionKind.PI_PI,
        InteractionKind.CATION_PI,
        InteractionKind.HOST_GUEST,
        InteractionKind.ION_PAIR,
    ]
    for idx, kind in enumerate(kinds):
        a = graph.add_atom("C", float(idx) * 100.0, 0.0).id
        b = graph.add_atom("C", float(idx) * 100.0 + 40.0, 0.0).id
        graph.add_bond(a, b, style=BondStyle.INTERACTION, interaction_kind=kind.value)

    layers = build_multilayer_chemical_graph(graph)

    assert {edge.kind for edge in layers.interaction_graph.edges} == set(kinds)
    assert sum(
        1
        for constraint in layers.layout_constraint_graph.constraints
        if constraint.kind == LayoutConstraintKind.INTERACTION_DISTANCE
    ) == len(kinds)
