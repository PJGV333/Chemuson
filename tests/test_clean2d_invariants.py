from __future__ import annotations

import math

import pytest


from chemuson.clean2d import (
    Clean2DInvariantError,
    assert_clean2d_invariants,
    capture_clean2d_snapshot,
)
from chemuson.core.model import BondStereo, BondStyle, MolGraph


def _metadata_graph() -> MolGraph:
    graph = MolGraph()
    c = graph.add_atom(
        "C",
        0.0,
        0.0,
        charge=-1,
        isotope=13,
        radical_electrons=1,
        explicit_h=1,
        stereo_cip="R",
    )
    o = graph.add_atom("O", 40.0, 0.0, charge=1, group_h_cap=1)
    graph.add_bond(
        c.id,
        o.id,
        order=2,
        style=BondStyle.WEDGE,
        stereo=BondStereo.UP,
        stereo_ez="E",
        is_aromatic=False,
        display_order=2,
    )
    return graph


def test_clean2d_invariants_allow_coordinate_only_change() -> None:
    graph = _metadata_graph()
    before_snapshot = capture_clean2d_snapshot(graph)
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    after = {1: (4.0, 2.0), 2: (44.0, 2.0)}

    assert_clean2d_invariants(before_snapshot, graph, before, after)


def test_clean2d_invariants_reject_atom_metadata_change() -> None:
    before_graph = _metadata_graph()
    after_graph = _metadata_graph()
    after_graph.atoms[1].charge = 0
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in before_graph.atoms.items()}
    after = dict(before)

    with pytest.raises(Clean2DInvariantError, match="atomos"):
        assert_clean2d_invariants(before_graph, after_graph, before, after)


def test_clean2d_invariants_reject_bond_metadata_change() -> None:
    before_graph = _metadata_graph()
    after_graph = _metadata_graph()
    after_graph.bonds[1].is_aromatic = True
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in before_graph.atoms.items()}
    after = dict(before)

    with pytest.raises(Clean2DInvariantError, match="enlaces"):
        assert_clean2d_invariants(before_graph, after_graph, before, after)


def test_clean2d_invariants_reject_missing_and_nan_coords() -> None:
    graph = _metadata_graph()
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}

    with pytest.raises(Clean2DInvariantError, match="faltante"):
        assert_clean2d_invariants(graph, graph, before, {1: (0.0, 0.0)})

    with pytest.raises(Clean2DInvariantError, match="no_finita"):
        assert_clean2d_invariants(graph, graph, before, {1: (0.0, 0.0), 2: (math.nan, 0.0)})
