"""Pruebas del motor puro clean2d_v2."""

from __future__ import annotations

import math
import os
import sys

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


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _center(points: list[tuple[float, float]]) -> tuple[float, float]:
    return (
        sum(x for x, _y in points) / len(points),
        sum(y for _x, y in points) / len(points),
    )
