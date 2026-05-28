"""Pruebas del motor puro clean2d_v2."""

from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import Clean2DParameters, optimize_clean2d_positions
from chemuson.core.model import MolGraph


def test_clean2d_v2_shortens_stretched_bond_without_mutating_graph() -> None:
    graph = MolGraph()
    a = graph.add_atom("C", 0.0, 0.0)
    b = graph.add_atom("C", 160.0, 0.0)
    graph.add_bond(a.id, b.id, order=1)

    after = optimize_clean2d_positions(
        graph,
        params=Clean2DParameters.quick(target_bond_length=42.0),
    )

    assert graph.get_atom(b.id).x == 160.0
    before_len = 160.0
    after_len = _distance(after[a.id], after[b.id])
    assert after_len < before_len
    assert abs(after_len - 42.0) < abs(before_len - 42.0)


def test_clean2d_v2_separates_overlapping_nonbonded_atoms() -> None:
    graph = MolGraph()
    a = graph.add_atom("C", 0.0, 0.0)
    b = graph.add_atom("C", 42.0, 0.0)
    c = graph.add_atom("O", 21.0, 0.0)
    graph.add_bond(a.id, b.id, order=1)

    after = optimize_clean2d_positions(
        graph,
        atom_ids={a.id, b.id, c.id},
        params=Clean2DParameters.publication(target_bond_length=42.0),
    )

    assert _distance(after[c.id], after[a.id]) > 21.0
    assert _distance(after[c.id], after[b.id]) > 21.0


def test_clean2d_v2_keeps_centroid_stable() -> None:
    graph = MolGraph()
    atoms = [
        graph.add_atom("C", 10.0, 10.0),
        graph.add_atom("C", 70.0, 20.0),
        graph.add_atom("C", 120.0, 80.0),
    ]
    graph.add_bond(atoms[0].id, atoms[1].id, order=1)
    graph.add_bond(atoms[1].id, atoms[2].id, order=1)

    before_center = _center([(atom.x, atom.y) for atom in atoms])
    after = optimize_clean2d_positions(
        graph,
        params=Clean2DParameters.quick(target_bond_length=42.0),
    )
    after_center = _center(list(after.values()))

    assert math.isclose(before_center[0], after_center[0], abs_tol=1e-6)
    assert math.isclose(before_center[1], after_center[1], abs_tol=1e-6)


def test_clean2d_v2_preserves_existing_good_bond_angles() -> None:
    graph = MolGraph()
    center = graph.add_atom("C", 0.0, 0.0)
    left = graph.add_atom("C", -36.37, -21.0)
    right = graph.add_atom("C", 36.37, -21.0)
    up = graph.add_atom("C", 0.0, 42.0)
    graph.add_bond(center.id, left.id, order=1)
    graph.add_bond(center.id, right.id, order=1)
    graph.add_bond(center.id, up.id, order=1)

    before_angles = sorted(
        _angle_between(
            (graph.get_atom(center.id).x, graph.get_atom(center.id).y),
            (graph.get_atom(a.id).x, graph.get_atom(a.id).y),
            (graph.get_atom(b.id).x, graph.get_atom(b.id).y),
        )
        for a, b in ((left, right), (left, up), (right, up))
    )
    after = optimize_clean2d_positions(
        graph,
        params=Clean2DParameters.quick(target_bond_length=42.0),
    )
    after_angles = sorted(
        _angle_between(after[center.id], after[a.id], after[b.id])
        for a, b in ((left, right), (left, up), (right, up))
    )

    assert after_angles == pytest.approx(before_angles, abs=1.0)


def test_clean2d_v2_does_not_create_large_local_deformations_on_chain() -> None:
    graph = MolGraph()
    atoms = [
        graph.add_atom("C", 0.0, 0.0),
        graph.add_atom("C", 42.0, 0.0),
        graph.add_atom("C", 63.0, 36.37),
        graph.add_atom("C", 105.0, 36.37),
        graph.add_atom("O", 126.0, 72.74),
    ]
    for left, right in zip(atoms, atoms[1:]):
        graph.add_bond(left.id, right.id, order=1)

    after = optimize_clean2d_positions(
        graph,
        params=Clean2DParameters.publication(target_bond_length=42.0),
    )

    max_move = max(
        _distance((atom.x, atom.y), after[atom.id])
        for atom in atoms
    )
    assert max_move < 8.0


def test_clean2d_v2_equalizes_distorted_acyclic_chain_lengths_and_angles() -> None:
    graph = MolGraph()
    atoms = [
        graph.add_atom("C", 0.0, 0.0),
        graph.add_atom("C", 80.0, 0.0),
        graph.add_atom("C", 120.0, 0.0),
        graph.add_atom("C", 200.0, 0.0),
    ]
    for left, right in zip(atoms, atoms[1:]):
        graph.add_bond(left.id, right.id, order=1)

    after = optimize_clean2d_positions(
        graph,
        params=Clean2DParameters.publication(target_bond_length=40.0),
    )

    lengths = [
        _distance(after[left.id], after[right.id])
        for left, right in zip(atoms, atoms[1:])
    ]
    assert lengths == pytest.approx([40.0, 40.0, 40.0], abs=1.0)

    internal_angles = [
        _angle_between(after[atoms[1].id], after[atoms[0].id], after[atoms[2].id]),
        _angle_between(after[atoms[2].id], after[atoms[1].id], after[atoms[3].id]),
    ]
    assert max(internal_angles) < 145.0


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _center(points: list[tuple[float, float]]) -> tuple[float, float]:
    return (
        sum(x for x, _y in points) / len(points),
        sum(y for _x, y in points) / len(points),
    )


def _angle_between(
    center: tuple[float, float],
    left: tuple[float, float],
    right: tuple[float, float],
) -> float:
    a1 = math.atan2(left[1] - center[1], left[0] - center[0])
    a2 = math.atan2(right[1] - center[1], right[0] - center[0])
    diff = abs((math.degrees(a2 - a1) + 180.0) % 360.0 - 180.0)
    return diff
