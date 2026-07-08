from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from chemuson.clean2d import Clean2DMode
from chemuson.core.model import BondStereo, BondStyle, MolGraph


Target = tuple[int, ...] | None


@dataclass(frozen=True)
class Clean2DRegressionCase:
    name: str
    family: str
    builder: Callable[[], MolGraph]
    mode: Clean2DMode = Clean2DMode.QUICK
    target: Target = None
    expected_states: frozenset[str] = frozenset({"applied", "no-op", "preserve-only"})
    tags: tuple[str, ...] = ("baseline",)
    known_delicate: bool = False
    known_failure: bool = False
    notes: str = ""

    @property
    def target_label(self) -> str:
        return "whole" if self.target is None else "selection"

    def has_tag(self, tag: str) -> bool:
        return tag in self.tags


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


def _substituted_benzene_graph(substituents: dict[int, tuple[str, float, float]]) -> MolGraph:
    graph = _benzene_graph()
    next_atom = 7
    next_bond = 7
    for ring_atom, (element, dx, dy) in substituents.items():
        atom = graph.add_atom(element, graph.atoms[ring_atom].x + dx, graph.atoms[ring_atom].y + dy, atom_id=next_atom)
        graph.add_bond(ring_atom, atom.id, bond_id=next_bond)
        next_atom += 1
        next_bond += 1
    return graph


def _toluene_like_graph() -> MolGraph:
    return _substituted_benzene_graph({1: ("C", 38.0, 22.0)})


def _phenol_like_graph() -> MolGraph:
    return _substituted_benzene_graph({1: ("O", 38.0, 22.0)})


def _aniline_like_graph() -> MolGraph:
    return _substituted_benzene_graph({1: ("N", 38.0, 22.0)})


def _para_disubstituted_benzene_graph() -> MolGraph:
    return _substituted_benzene_graph({1: ("O", 38.0, 22.0), 4: ("C", -38.0, -22.0)})


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


def _anthracene_like_graph() -> MolGraph:
    graph = _fused_ring_graph()
    for atom_id, (x, y) in {
        11: (180.0, 0.0),
        12: (216.0, 24.0),
        13: (216.0, 66.0),
        14: (180.0, 90.0),
    }.items():
        graph.add_atom("C", x, y, atom_id=atom_id)
    for bond_id, (a1, a2) in enumerate(((8, 11), (11, 12), (12, 13), (13, 14), (14, 9)), start=12):
        graph.add_bond(a1, a2, bond_id=bond_id, order=1, is_aromatic=True)
    return graph


def _pyridine_like_graph() -> MolGraph:
    graph = _benzene_graph()
    graph.atoms[1].element = "N"
    return graph


def _five_member_heteroaromatic_graph(elements: tuple[str, ...]) -> MolGraph:
    graph = MolGraph()
    radius = 36.0
    for idx, element in enumerate(elements):
        angle = math.radians(idx * 72.0 - 90.0)
        graph.add_atom(element, math.cos(angle) * radius, math.sin(angle) * radius, atom_id=idx + 1)
    for idx in range(len(elements)):
        graph.add_bond(idx + 1, ((idx + 1) % len(elements)) + 1, bond_id=idx + 1, order=1, is_aromatic=True)
    return graph


def _imidazole_like_graph() -> MolGraph:
    return _five_member_heteroaromatic_graph(("N", "C", "N", "C", "C"))


def _furan_like_graph() -> MolGraph:
    return _five_member_heteroaromatic_graph(("O", "C", "C", "C", "C"))


def _thiophene_like_graph() -> MolGraph:
    return _five_member_heteroaromatic_graph(("S", "C", "C", "C", "C"))


def _ammonium_like_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("N", 0.0, 0.0, atom_id=1, charge=1, explicit_h=4)
    graph.add_atom("C", 42.0, 0.0, atom_id=2)
    graph.add_bond(1, 2, bond_id=1)
    return graph


def _carboxylate_like_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    graph.add_atom("O", 42.0, 24.0, atom_id=2)
    graph.add_atom("O", 42.0, -24.0, atom_id=3, charge=-1)
    graph.add_atom("C", -42.0, 0.0, atom_id=4)
    graph.add_bond(1, 2, bond_id=1, order=2)
    graph.add_bond(1, 3, bond_id=2, order=1)
    graph.add_bond(1, 4, bond_id=3, order=1)
    return graph


def _pyridinium_like_graph() -> MolGraph:
    graph = _pyridine_like_graph()
    graph.atoms[1].charge = 1
    return graph


def _zwitterion_like_graph() -> MolGraph:
    graph = _carboxylate_like_graph()
    graph.add_atom("N", -84.0, 0.0, atom_id=5, charge=1, explicit_h=3)
    graph.add_bond(4, 5, bond_id=4)
    return graph


def _chiral_wedge_hash_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1, stereo_cip="R")
    graph.add_atom("F", 42.0, 0.0, atom_id=2)
    graph.add_atom("Cl", -42.0, 0.0, atom_id=3)
    graph.add_atom("Br", 0.0, 42.0, atom_id=4)
    graph.add_atom("H", 0.0, -42.0, atom_id=5, is_explicit=True)
    graph.add_bond(1, 2, bond_id=1, style=BondStyle.WEDGE, stereo=BondStereo.UP)
    graph.add_bond(1, 3, bond_id=2, style=BondStyle.HASHED, stereo=BondStereo.DOWN)
    graph.add_bond(1, 4, bond_id=3)
    graph.add_bond(1, 5, bond_id=4)
    return graph


def _alkene_ez_like_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0, atom_id=1)
    graph.add_atom("C", 42.0, 0.0, atom_id=2)
    graph.add_atom("Cl", -24.0, 30.0, atom_id=3)
    graph.add_atom("H", -24.0, -30.0, atom_id=4, is_explicit=True)
    graph.add_atom("Br", 66.0, -30.0, atom_id=5)
    graph.add_atom("H", 66.0, 30.0, atom_id=6, is_explicit=True)
    graph.add_bond(1, 2, bond_id=1, order=2, stereo_ez="E")
    graph.add_bond(1, 3, bond_id=2)
    graph.add_bond(1, 4, bond_id=3)
    graph.add_bond(2, 5, bond_id=4)
    graph.add_bond(2, 6, bond_id=5)
    return graph


def _partial_selection_long_chain_graph() -> MolGraph:
    graph = MolGraph()
    for idx in range(6):
        graph.add_atom("C", idx * 42.0, float((idx % 2) * 4), atom_id=idx + 1)
    for idx in range(5):
        graph.add_bond(idx + 1, idx + 2, bond_id=idx + 1)
    return graph


def _partial_selection_ring_touch_graph() -> MolGraph:
    graph = _toluene_like_graph()
    graph.add_atom("O", graph.atoms[7].x + 42.0, graph.atoms[7].y, atom_id=8)
    graph.add_bond(7, 8, bond_id=8)
    return graph


def _biphenyl_like_graph() -> MolGraph:
    graph = _benzene_graph()
    offset = 118.0
    for idx in range(6):
        angle = math.radians(idx * 60.0 + 30.0)
        graph.add_atom("C", offset + math.cos(angle) * 42.0, math.sin(angle) * 42.0, atom_id=idx + 7)
    for idx in range(6):
        graph.add_bond(idx + 7, ((idx + 1) % 6) + 7, bond_id=idx + 7, order=1, is_aromatic=True)
    graph.add_bond(1, 10, bond_id=13)
    return graph


def _diphenyl_ether_like_graph() -> MolGraph:
    graph = _biphenyl_like_graph()
    graph.remove_bond(13)
    graph.add_atom("O", 59.0, 20.0, atom_id=13)
    graph.add_bond(1, 13, bond_id=13)
    graph.add_bond(13, 10, bond_id=14)
    return graph


def _triphenyl_like_graph() -> MolGraph:
    graph = _biphenyl_like_graph()
    for idx in range(6):
        angle = math.radians(idx * 60.0 + 30.0)
        graph.add_atom("C", 59.0 + math.cos(angle) * 42.0, 108.0 + math.sin(angle) * 42.0, atom_id=idx + 13)
    for idx in range(6):
        graph.add_bond(idx + 13, ((idx + 1) % 6) + 13, bond_id=idx + 14, order=1, is_aromatic=True)
    graph.add_bond(4, 16, bond_id=20)
    return graph


def _macrocycle_like_graph() -> MolGraph:
    graph = MolGraph()
    radius = 82.0
    for idx in range(12):
        angle = math.radians(idx * 30.0)
        graph.add_atom("C", math.cos(angle) * radius, math.sin(angle) * radius, atom_id=idx + 1)
    for idx in range(12):
        graph.add_bond(idx + 1, ((idx + 1) % 12) + 1, bond_id=idx + 1)
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
        tags=("baseline",),
    ),
    Clean2DRegressionCase(
        name="aromatic_benzene_regular",
        family="aromatic",
        builder=_benzene_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "known_delicate"),
        known_delicate=True,
        notes="Current aromatic baseline may reject candidates as worse-quality; no geometry fix in this change.",
    ),
    Clean2DRegressionCase(
        name="fused_aromatic_current_baseline",
        family="fused-ring",
        builder=_fused_ring_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
        notes="Observational fused-ring baseline; this change does not improve fused geometry.",
    ),
    Clean2DRegressionCase(
        name="partial_selection_boundary_chain",
        family="partial-selection",
        builder=_partial_selection_graph,
        target=(2, 3),
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("known_delicate", "selection_boundary", "stereo_sensitive"),
        known_delicate=True,
        notes="Records boundary-sensitive current behaviour without changing boundary policy.",
    ),
    Clean2DRegressionCase(
        name="dirty_coordinate_valid_graph",
        family="dirty-import",
        builder=_dirty_coordinate_graph,
        expected_states=frozenset({"applied", "failed-controlled"}),
        tags=("known_delicate",),
        known_delicate=True,
        notes="Valid graph with poor starting coordinates; geometry quality is not fixed here.",
    ),
    Clean2DRegressionCase("aromatic_toluene_like", "substituted-aromatic", _toluene_like_graph, tags=("baseline",)),
    Clean2DRegressionCase("aromatic_phenol_like", "substituted-aromatic", _phenol_like_graph, tags=("baseline",)),
    Clean2DRegressionCase("aromatic_aniline_like", "substituted-aromatic", _aniline_like_graph, tags=("baseline",)),
    Clean2DRegressionCase("aromatic_para_disubstituted_like", "substituted-aromatic", _para_disubstituted_benzene_graph, tags=("baseline",)),
    Clean2DRegressionCase(
        "fused_aromatic_anthracene_like",
        "fused-ring",
        _anthracene_like_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "heteroaromatic_pyridine_like",
        "heteroaromatic",
        _pyridine_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "heteroaromatic_imidazole_like",
        "heteroaromatic",
        _imidazole_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "heteroaromatic_furan_like",
        "heteroaromatic",
        _furan_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "heteroaromatic_thiophene_like",
        "heteroaromatic",
        _thiophene_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase("charged_ammonium_like", "charged", _ammonium_like_graph, tags=("baseline", "charged")),
    Clean2DRegressionCase(
        "charged_carboxylate_like",
        "charged",
        _carboxylate_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "charged", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "charged_pyridinium_like",
        "charged",
        _pyridinium_like_graph,
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "charged", "known_delicate"),
        known_delicate=True,
    ),
    Clean2DRegressionCase("charged_zwitterion_like", "charged", _zwitterion_like_graph, tags=("baseline", "charged")),
    Clean2DRegressionCase("stereo_chiral_wedge_hash", "stereo-sensitive", _chiral_wedge_hash_graph, tags=("baseline", "stereo_sensitive")),
    Clean2DRegressionCase("stereo_alkene_ez_like", "stereo-sensitive", _alkene_ez_like_graph, tags=("baseline", "stereo_sensitive")),
    Clean2DRegressionCase(
        "selection_boundary_acyclic_internal",
        "selection-boundary",
        _partial_selection_long_chain_graph,
        target=(3, 4),
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("baseline", "selection_boundary"),
    ),
    Clean2DRegressionCase(
        "selection_boundary_ring_substituent",
        "selection-boundary",
        _partial_selection_ring_touch_graph,
        target=(1, 7),
        expected_states=frozenset({"applied", "no-op", "failed-controlled"}),
        tags=("known_delicate", "selection_boundary"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "selection_boundary_biphenyl_linker",
        "selection-boundary",
        _biphenyl_like_graph,
        target=(1, 10),
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "selection_boundary", "complex_policy_guard"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "multiblock_biphenyl_like",
        "multi-block",
        _biphenyl_like_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "multiblock_diphenyl_ether_like",
        "multi-block",
        _diphenyl_ether_like_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "multiblock_triphenyl_like",
        "multi-block",
        _triphenyl_like_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
    ),
    Clean2DRegressionCase(
        "macrocycle_twelve_member_baseline",
        "macrocycle",
        _macrocycle_like_graph,
        expected_states=frozenset({"applied", "no-op", "preserve-only", "failed-controlled"}),
        tags=("known_delicate", "complex_policy_guard"),
        known_delicate=True,
    ),
)


def get_regression_cases() -> tuple[Clean2DRegressionCase, ...]:
    return REGRESSION_CASES
