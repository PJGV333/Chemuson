from __future__ import annotations

import math
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.clean2d import run_clean2d_engine
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
