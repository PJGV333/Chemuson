from __future__ import annotations

import math

import pytest


from chemuson.clean2d import (
    Clean2DCandidate,
    clean2d_geometry_hash,
    generate_clean2d_candidates,
    rank_clean2d_candidates,
    ring_degeneracy_score,
    run_clean2d_engine,
)
from chemuson.chemio.rdkit_safe import _project_missing_clean2d_hydrogens
from chemuson.core.model import MolGraph


def _regular_ring(
    graph: MolGraph,
    *,
    n: int,
    center: tuple[float, float] = (0.0, 0.0),
    bond_length: float = 40.0,
    aromatic: bool = False,
    element: str = "C",
) -> list[int]:
    radius = bond_length / (2.0 * math.sin(math.pi / n))
    ids: list[int] = []
    for idx in range(n):
        angle = math.radians(360.0 * idx / n + (30.0 if n == 6 else -90.0))
        atom = graph.add_atom(
            element,
            center[0] + radius * math.cos(angle),
            center[1] + radius * math.sin(angle),
        )
        ids.append(atom.id)
    for idx in range(n):
        graph.add_bond(ids[idx], ids[(idx + 1) % n], order=1, is_aromatic=aromatic)
    return ids


def _distance(graph: MolGraph, a_id: int, b_id: int, coords: dict[int, tuple[float, float]]) -> float:
    return math.hypot(coords[b_id][0] - coords[a_id][0], coords[b_id][1] - coords[a_id][1])


def _angle(
    coords: dict[int, tuple[float, float]],
    left: int,
    center: int,
    right: int,
) -> float:
    a1 = math.atan2(coords[left][1] - coords[center][1], coords[left][0] - coords[center][0])
    a2 = math.atan2(coords[right][1] - coords[center][1], coords[right][0] - coords[center][0])
    return abs((math.degrees(a2 - a1) + 180.0) % 360.0 - 180.0)


def test_generate_candidates_attempts_rdkit_for_cyclic_graphs() -> None:
    graph = MolGraph()
    ring = []
    for idx, (x, y) in enumerate([(0, 0), (4, 0), (5, 1), (3, 2), (1, 2), (-1, 1)], 1):
        ring.append(graph.add_atom("C", float(x), float(y), atom_id=idx).id)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)
    candidates = generate_clean2d_candidates(graph, ring, mode="publication", target_bond_length=40.0)
    sources = {candidate.source for candidate in candidates}

    assert "rdkit_isolated" in sources
    assert "internal_templates" in sources


def test_engine_rebuilds_distorted_ring_with_internal_candidate() -> None:
    graph = MolGraph()
    ring = []
    for idx, (x, y) in enumerate([(0, 0), (4, 0), (5, 1), (3, 2), (1, 2), (-1, 1)], 1):
        ring.append(graph.add_atom("C", float(x), float(y), atom_id=idx).id)
    for idx in range(6):
        graph.add_bond(ring[idx], ring[(idx + 1) % 6], order=1, is_aromatic=True)

    result = run_clean2d_engine(graph, ring, mode="publication", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert result.selected.source != "current"
    assert ring_degeneracy_score(result.selected.coords, set(ring)) > 0.25
    lengths = [
        _distance(graph, ring[idx], ring[(idx + 1) % 6], result.selected.coords)
        for idx in range(6)
    ]
    assert sum(lengths) / len(lengths) == pytest.approx(40.0, rel=0.08)


def test_quick_mode_keeps_already_clean_ring_stable() -> None:
    graph = MolGraph()
    ring = _regular_ring(graph, n=6, aromatic=True)

    result = run_clean2d_engine(graph, ring, mode="quick", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    assert result.selected.source == "current"
    assert result.selected.novelty < 0.5


def test_engine_rebuilds_linear_sp3_chain_as_zigzag() -> None:
    graph = MolGraph()
    atoms = [graph.add_atom("C", float(idx * 40), 0.0) for idx in range(4)]
    for left, right in zip(atoms, atoms[1:]):
        graph.add_bond(left.id, right.id, order=1)

    result = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)

    assert result.ok
    assert result.selected is not None
    coords = result.selected.coords
    assert _angle(coords, atoms[0].id, atoms[1].id, atoms[2].id) < 155.0
    assert _angle(coords, atoms[1].id, atoms[2].id, atoms[3].id) < 155.0


def test_rdkit_safe_projection_completes_explicit_hydrogens() -> None:
    graph = MolGraph()
    c = graph.add_atom("C", 0.0, 0.0)
    h1 = graph.add_atom("H", 1.0, 0.0)
    h2 = graph.add_atom("H", 0.0, 1.0)
    graph.add_bond(c.id, h1.id, order=1)
    graph.add_bond(c.id, h2.id, order=1)

    projected = _project_missing_clean2d_hydrogens(graph, {c.id: (10.0, 20.0)})

    assert projected[h1.id] == pytest.approx((11.0, 20.0), rel=1e-9)
    assert projected[h2.id] == pytest.approx((10.0, 21.0), rel=1e-9)


def test_avoid_hashes_only_blocks_clean_propose_alternatives() -> None:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0)
    a2 = graph.add_atom("C", 100.0, 0.0)
    graph.add_bond(a1.id, a2.id, order=1)
    repaired = {a1.id: (0.0, 0.0), a2.id: (40.0, 0.0)}
    repaired_hash = clean2d_geometry_hash(graph, repaired)
    candidate = Clean2DCandidate(
        source="internal_templates",
        coords=repaired,
        message="repair",
    )
    before_bad = {a1.id: (0.0, 0.0), a2.id: (100.0, 0.0)}

    quick = rank_clean2d_candidates(
        graph,
        [candidate],
        before_bad,
        {a1.id, a2.id},
        mode="quick",
        target_bond_length=40.0,
        avoid_hashes={repaired_hash},
    )
    publication = rank_clean2d_candidates(
        graph,
        [candidate],
        before_bad,
        {a1.id, a2.id},
        mode="publication",
        target_bond_length=40.0,
        avoid_hashes={repaired_hash},
    )
    propose_bad = rank_clean2d_candidates(
        graph,
        [candidate],
        before_bad,
        {a1.id, a2.id},
        mode="propose",
        target_bond_length=40.0,
        avoid_hashes={repaired_hash},
    )
    before_clean = dict(repaired)
    propose_clean = rank_clean2d_candidates(
        graph,
        [candidate],
        before_clean,
        {a1.id, a2.id},
        mode="propose",
        target_bond_length=40.0,
        avoid_hashes={repaired_hash},
    )

    assert quick.ok
    assert publication.ok
    assert not propose_bad.ok
    assert propose_bad.message == "La estructura debe optimizarse antes de proponer conformeros 2D"
    assert propose_clean.selected is None
    assert any(item.rejection_reason == "geometria_repetida" for item in propose_clean.rejected)
