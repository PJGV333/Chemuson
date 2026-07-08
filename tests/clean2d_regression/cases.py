from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from chemuson.clean2d import Clean2DMode
from chemuson.core.model import BondStyle, MolGraph


Target = tuple[int, ...] | None


@dataclass(frozen=True)
class Clean2DRegressionCase:
    name: str
    family: str
    builder: Callable[[], MolGraph]
    mode: Clean2DMode = Clean2DMode.QUICK
    target: Target = None
    expected_states: frozenset[str] = frozenset({"applied", "no-op", "preserve-only"})
    known_delicate: bool = False
    known_failure: bool = False
    notes: str = ""

    @property
    def target_label(self) -> str:
        return "whole" if self.target is None else "selection"


def _chain_graph() -> MolGraph:
    graph = MolGraph()
    for idx, (x, y) in enumerate(((0.0, 0.0), (64.0, 3.0), (132.0, -4.0), (198.0, 2.0)), start=1):
        graph.add_atom("C", x, y, atom_id=idx)
    graph.add_bond(1, 2, bond_id=1)
    graph.add_bond(2, 3, bond_id=2)
    graph.add_bond(3, 4, bond_id=3)
    return graph


def _benzene_graph() -> MolGraph:
    graph = MolGraph()
    radius = 42.0
    for idx in range(6):
        angle = math.radians(idx * 60.0 + 30.0)
        graph.add_atom("C", math.cos(angle) * radius, math.sin(angle) * radius, atom_id=idx + 1)
    for idx in range(6):
        graph.add_bond(idx + 1, ((idx + 1) % 6) + 1, bond_id=idx + 1, order=1, is_aromatic=True)
    return graph


def _fused_ring_graph() -> MolGraph:
    graph = MolGraph()
    coords = {
        1: (0.0, 24.0),
        2: (36.0, 0.0),
        3: (72.0, 24.0),
        4: (72.0, 66.0),
        5: (36.0, 90.0),
        6: (0.0, 66.0),
        7: (108.0, 0.0),
        8: (144.0, 24.0),
        9: (144.0, 66.0),
        10: (108.0, 90.0),
    }
    for atom_id, (x, y) in coords.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    for bond_id, (a1, a2) in enumerate(
        ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1), (3, 7), (7, 8), (8, 9), (9, 10), (10, 4)),
        start=1,
    ):
        graph.add_bond(a1, a2, bond_id=bond_id, order=1, is_aromatic=True)
    return graph


def _partial_selection_graph() -> MolGraph:
    graph = _chain_graph()
    graph.atoms[2].stereo_cip = "R"
    graph.bonds[2].style = BondStyle.WEDGE
    return graph


def _dirty_coordinate_graph() -> MolGraph:
    graph = MolGraph()
    for atom_id, (x, y) in {
        1: (0.0, 0.0),
        2: (9.0, 0.5),
        3: (18.0, 0.0),
        4: (9.0, 1.0),
        5: (26.0, 0.4),
    }.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    graph.add_bond(1, 2, bond_id=1)
    graph.add_bond(2, 3, bond_id=2, order=2, stereo_ez="E")
    graph.add_bond(3, 4, bond_id=3)
    graph.add_bond(4, 5, bond_id=4)
    return graph


REGRESSION_CASES: tuple[Clean2DRegressionCase, ...] = (
    Clean2DRegressionCase(
        name="acyclic_butane_stretched",
        family="acyclic",
        builder=_chain_graph,
        expected_states=frozenset({"applied", "no-op"}),
    ),
    Clean2DRegressionCase(
        name="aromatic_benzene_regular",
        family="aromatic",
        builder=_benzene_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        known_delicate=True,
        notes="Current aromatic baseline may reject candidates as worse-quality; no geometry fix in this change.",
    ),
    Clean2DRegressionCase(
        name="fused_aromatic_current_baseline",
        family="fused-ring",
        builder=_fused_ring_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        known_delicate=True,
        notes="Observational fused-ring baseline; this change does not improve fused geometry.",
    ),
    Clean2DRegressionCase(
        name="partial_selection_boundary_chain",
        family="partial-selection",
        builder=_partial_selection_graph,
        target=(2, 3),
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        known_delicate=True,
        notes="Records boundary-sensitive current behaviour without changing boundary policy.",
    ),
    Clean2DRegressionCase(
        name="dirty_coordinate_valid_graph",
        family="dirty-import",
        builder=_dirty_coordinate_graph,
        expected_states=frozenset({"applied", "failed-controlled"}),
        known_delicate=True,
        notes="Valid graph with poor starting coordinates; geometry quality is not fixed here.",
    ),
)


def get_regression_cases() -> tuple[Clean2DRegressionCase, ...]:
    return REGRESSION_CASES
