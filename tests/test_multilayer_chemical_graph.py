from __future__ import annotations



from chemuson.core.layers import (
    BlockEdgeKind,
    BlockKind,
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


def test_block_graph_detects_fused_rings_linkers_substituents_and_stereo() -> None:
    graph = MolGraph()
    atoms = {idx: graph.add_atom("C", float(idx) * 10.0, 0.0).id for idx in range(1, 24)}
    for left, right in [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]:
        graph.add_bond(atoms[left], atoms[right], is_aromatic=True)
    for left, right in [(5, 10), (10, 9), (9, 8), (8, 7), (7, 6)]:
        graph.add_bond(atoms[left], atoms[right], is_aromatic=True)
    for left, right in [(17, 18), (18, 19), (19, 20), (20, 21), (21, 22), (22, 17)]:
        graph.add_bond(atoms[left], atoms[right], is_aromatic=True)
    graph.add_bond(atoms[3], atoms[11], order=1)
    graph.add_bond(atoms[11], atoms[12], order=1)
    graph.add_bond(atoms[12], atoms[17], order=1)
    graph.add_bond(atoms[8], atoms[13], order=1)
    graph.add_atom("C", 200.0, 0.0, atom_id=16, stereo_cip="R")
    graph.add_bond(atoms[2], 16, order=1)

    layers = build_multilayer_chemical_graph(graph)
    kinds = [block.kind for block in layers.block_graph.blocks]

    assert BlockKind.FUSED_SYSTEM in kinds
    assert BlockKind.AROMATIC_RING in kinds
    assert BlockKind.LINKER in kinds
    assert BlockKind.TERMINAL_SUBSTITUENT in kinds
    assert BlockKind.STEREO_CENTER in kinds
    assert any(edge.kind == BlockEdgeKind.LINKER for edge in layers.block_graph.edges)


def test_block_graph_detects_macrocycle_bridge_and_internal_cavity() -> None:
    graph = MolGraph()
    atoms = [graph.add_atom("C", float(idx), 0.0).id for idx in range(12)]
    for idx, atom_id in enumerate(atoms):
        graph.add_bond(atom_id, atoms[(idx + 1) % len(atoms)], order=1)
    bridge_a = graph.add_atom("C", 30.0, 20.0).id
    bridge_b = graph.add_atom("C", 40.0, 20.0).id
    graph.add_bond(atoms[2], bridge_a, order=1)
    graph.add_bond(bridge_a, bridge_b, order=1)
    graph.add_bond(bridge_b, atoms[8], order=1)

    layers = build_multilayer_chemical_graph(graph)
    kinds = [block.kind for block in layers.block_graph.blocks]

    assert BlockKind.MACROCYCLE in kinds
    assert BlockKind.INTRAMOLECULAR_BRIDGE in kinds
    assert BlockKind.INTERNAL_CAVITY in kinds
