from __future__ import annotations

import math


from chemuson.clean2d import clean2d_geometry_hash, run_clean2d_engine
from chemuson.core.model import MolGraph


def _add_benzene(graph: MolGraph, cx: float, cy: float) -> list[int]:
    ids: list[int] = []
    radius = 40.0
    for idx in range(6):
        angle = math.radians(60.0 * idx + 30.0)
        ids.append(graph.add_atom("C", cx + radius * math.cos(angle), cy + radius * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(ids[idx], ids[(idx + 1) % 6], order=1, is_aromatic=True)
    return ids


def _biphenyl() -> MolGraph:
    graph = MolGraph()
    left = _add_benzene(graph, 0.0, 0.0)
    right = _add_benzene(graph, 109.3, 40.0)
    graph.add_bond(left[0], right[3], order=1)
    return graph


def test_quick_keeps_clean_biphenyl_but_propose_returns_alternative() -> None:
    graph = _biphenyl()
    before = {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}
    before_hash = clean2d_geometry_hash(graph, before)

    quick = run_clean2d_engine(graph, mode="quick", target_bond_length=40.0)
    propose = run_clean2d_engine(graph, mode="propose", target_bond_length=40.0)

    assert quick.ok
    assert quick.selected is not None
    assert quick.selected.source == "current"
    assert propose.ok
    assert propose.selected is not None
    assert propose.selected.message == "Estructura 2D alternativa propuesta"
    assert propose.selected.geometry_hash != before_hash
    assert propose.selected.novelty > 2.0


def test_propose_avoids_recent_geometry_hash_when_alternatives_exist() -> None:
    graph = _biphenyl()
    first = run_clean2d_engine(graph, mode="propose", target_bond_length=40.0)
    assert first.ok
    assert first.selected is not None

    second = run_clean2d_engine(
        graph,
        mode="propose",
        target_bond_length=40.0,
        avoid_hashes={first.selected.geometry_hash},
    )

    assert second.ok
    assert second.selected is not None
    assert second.selected.geometry_hash != first.selected.geometry_hash
