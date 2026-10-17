"""Reproducible molecular/geometry cases for the depiction visual lab."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import math
from typing import Callable

from chemuson.core.model import BondStereo, BondStyle, MolGraph
from chemuson.gui.templates import build_haworth_template


@dataclass(frozen=True)
class DepictionCase:
    case_id: str
    title: str
    category: str
    kind: str
    description: str
    build: Callable[[], MolGraph]
    editorial_calibration_eligible: bool = False
    editorial_exclusion_reason: str | None = None


def _atom(
    graph: MolGraph,
    element: str,
    x: float,
    y: float,
    *,
    explicit: bool | None = None,
) -> int:
    item = graph.add_atom(
        element,
        float(x),
        float(y),
        is_explicit=(element != "C" if explicit is None else explicit),
    )
    return item.id


def _bond(
    graph: MolGraph,
    a1: int,
    a2: int,
    *,
    order: int = 1,
    style: BondStyle = BondStyle.PLAIN,
    stereo: BondStereo = BondStereo.NONE,
    aromatic: bool = False,
    ring_id: int | None = None,
) -> None:
    graph.add_bond(
        a1,
        a2,
        order=order,
        style=style,
        stereo=stereo,
        is_aromatic=aromatic,
        ring_id=ring_id,
    )


def _chiral_label_wedge(label: str, wedge_style: BondStyle) -> MolGraph:
    graph = MolGraph()
    center = _atom(graph, "C", 0, 0)
    terminal = _atom(graph, label, 78, 0)
    _bond(
        graph,
        center,
        terminal,
        style=wedge_style,
        stereo=BondStereo.UP if wedge_style == BondStyle.WEDGE else BondStereo.DOWN,
    )
    branches = (
        (("OH", -42, -66), ("Cl", -42, 66), ("C", 0, 88))
        if label == "F"
        else (("F", -42, -66), ("Cl", -42, 66), ("C", 0, 88))
    )
    for element, x, y in branches:
        branch = _atom(graph, element, x, y)
        _bond(graph, center, branch)
    return graph


def _tetrahedral_pair() -> MolGraph:
    graph = MolGraph()
    center = _atom(graph, "C", 0, 0)
    groups = (
        ("F", -76, -10, BondStyle.WEDGE, BondStereo.UP),
        ("Cl", 68, 10, BondStyle.HASHED, BondStereo.DOWN),
        ("Br", 0, -82, BondStyle.PLAIN, BondStereo.NONE),
        ("C", 0, 86, BondStyle.PLAIN, BondStereo.NONE),
    )
    for element, x, y, style, stereo in groups:
        terminal = _atom(graph, element, x, y)
        _bond(graph, center, terminal, style=style, stereo=stereo)
    return graph


def _stereo_center_next_to_double(style: BondStyle) -> MolGraph:
    graph = MolGraph()
    center = _atom(graph, "C", 0, 0)
    fluorine = _atom(graph, "F", 78, 0)
    _bond(
        graph,
        center,
        fluorine,
        style=style,
        stereo=BondStereo.UP if style == BondStyle.WEDGE else BondStereo.DOWN,
    )
    alkene = _atom(graph, "C", -68, -40)
    _bond(graph, center, alkene, order=2)
    alkene_methyl = _atom(graph, "C", -142, -68)
    _bond(graph, alkene, alkene_methyl)
    alkene_chlorine = _atom(graph, "Cl", -142, -12)
    _bond(graph, alkene, alkene_chlorine)
    methyl = _atom(graph, "C", -68, 45)
    _bond(graph, center, methyl)
    return graph


def _carbonyl_oxygen() -> MolGraph:
    graph = MolGraph()
    carbonyl = _atom(graph, "C", 0, 0)
    oxygen = _atom(graph, "O", 82, 0)
    _bond(graph, carbonyl, oxygen, order=2)
    for y in (-48, 48):
        methyl = _atom(graph, "C", -78, y)
        _bond(graph, carbonyl, methyl)
    return graph


def _geometric_branch_fixture() -> MolGraph:
    """Abstract center for incompatible bond styles; deliberately not a molecule."""
    graph = MolGraph()
    center = _atom(graph, "C", 0, 0)
    branches = (
        ("F", 76, 0, BondStyle.WEDGE, 1, BondStereo.UP),
        ("C", -44, -68, BondStyle.BOLD, 2, BondStereo.NONE),
        ("C", -44, 68, BondStyle.PLAIN, 3, BondStereo.NONE),
    )
    for element, x, y, style, order, stereo in branches:
        end = _atom(graph, element, x, y)
        _bond(graph, center, end, order=order, style=style, stereo=stereo)
    return graph


def _ring(sides: int, bold_edges: tuple[int, ...]) -> MolGraph:
    graph = MolGraph()
    radius = 76.0
    atoms = [
        _atom(
            graph,
            "C",
            radius * math.cos(math.radians(-90 + index * 360.0 / sides)),
            radius * math.sin(math.radians(-90 + index * 360.0 / sides)),
        )
        for index in range(sides)
    ]
    for index, atom_id in enumerate(atoms):
        _bond(
            graph,
            atom_id,
            atoms[(index + 1) % sides],
            style=BondStyle.BOLD if index in bold_edges else BondStyle.PLAIN,
            ring_id=1,
        )
    return graph


def _haworth() -> MolGraph:
    return build_haworth_template(68.0, anomeric_up=True, bold_front=True)


def _aromatic_chiral_sidechain() -> MolGraph:
    graph = MolGraph()
    radius = 74.0
    ring = [
        _atom(
            graph,
            "C",
            radius * math.cos(math.radians(index * 60.0)),
            radius * math.sin(math.radians(index * 60.0)),
        )
        for index in range(6)
    ]
    for index, atom_id in enumerate(ring):
        order = 2 if index % 2 == 0 else 1
        _bond(
            graph,
            atom_id,
            ring[(index + 1) % 6],
            order=order,
            aromatic=True,
            ring_id=1,
        )
    attachment = ring[0]
    center = _atom(graph, "C", 116, -1)
    _bond(graph, attachment, center)
    fluorine = _atom(graph, "F", 178, -40)
    _bond(graph, center, fluorine, style=BondStyle.WEDGE, stereo=BondStereo.UP)
    oxygen = _atom(graph, "OH", 150, 64)
    _bond(graph, center, oxygen)
    methyl = _atom(graph, "C", 218, 15)
    _bond(graph, center, methyl)
    return graph


_EDITORIAL_EXCLUSION = "outside the approved A1–B3 editorial calibration cohort"

CASES: tuple[DepictionCase, ...] = (
    DepictionCase("A1", "Solid wedge · OH label", "A · Wedges and labels", "molecule", "Tetrahedral center with an OH-bearing wedge terminus.", lambda: _chiral_label_wedge("OH", BondStyle.WEDGE), True),
    DepictionCase("A2", "Solid wedge · F label", "A · Wedges and labels", "molecule", "Tetrahedral center with a compact halogen label.", lambda: _chiral_label_wedge("F", BondStyle.WEDGE), True),
    DepictionCase("A3", "Hashed wedge · OH label", "A · Wedges and labels", "molecule", "Tetrahedral center with a hashed OH-bearing substituent.", lambda: _chiral_label_wedge("OH", BondStyle.HASHED), True),
    DepictionCase("A4", "Tetrahedral center · wedge + hash", "A · Wedges and labels", "molecule", "Four distinct substituents; one solid and one hashed stereobond.", _tetrahedral_pair, True),
    DepictionCase("B1", "Hashed wedge + double bond", "B · Bond junctions", "molecule", "Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.", lambda: _stereo_center_next_to_double(BondStyle.HASHED), True),
    DepictionCase("B2", "Solid wedge + double bond", "B · Bond junctions", "molecule", "Tetrahedral stereocenter adjacent to a C=C sigma/pi depiction.", lambda: _stereo_center_next_to_double(BondStyle.WEDGE), True),
    DepictionCase("B3", "Carbonyl · oxygen label", "B · Bond junctions", "molecule", "Acetone carbonyl with the actual O label and neighboring sigma bonds.", _carbonyl_oxygen, True),
    DepictionCase("B4", "Three-way wedge / bold / triple", "B · Bond junctions", "geometric-fixture", "Geometric fixture / non-editorial: invalid combined bond-order/style stress case; excluded from preset selection and style calibration.", _geometric_branch_fixture, False, "geometric fixture; non-editorial; excluded from preset selection / style calibration"),
    DepictionCase("C1", "Cyclopentane · bold lower edges", "C · Rings and projections", "molecule", "Five-membered saturated ring with two visually bold edges.", lambda: _ring(5, (2, 3)), False, _EDITORIAL_EXCLUSION),
    DepictionCase("C2", "Cyclohexane · bold lower edges", "C · Rings and projections", "molecule", "Six-membered saturated ring with two visually bold edges.", lambda: _ring(6, (2, 3)), False, _EDITORIAL_EXCLUSION),
    DepictionCase("C3", "Haworth sugar template", "C · Rings and projections", "template-depiction", "Deterministic Haworth projection fixture; diagnostic only and not a source of AB3 calibration values.", _haworth, False, _EDITORIAL_EXCLUSION),
    DepictionCase("C4", "Aromatic ring · chiral side chain", "C · Rings and projections", "molecule", "Alternating aromatic ring with a tetrahedral fluorinated OH-bearing side chain.", _aromatic_chiral_sidechain, False, _EDITORIAL_EXCLUSION),
)

CASE_BY_ID = {case.case_id: case for case in CASES}


def case_graph(case_id: str) -> MolGraph:
    """Create a fresh graph for one stable gallery case ID."""
    try:
        case = CASE_BY_ID[case_id]
    except KeyError as exc:
        raise KeyError(f"unknown depiction case ID: {case_id}") from exc
    return case.build()


def clone_graph(graph: MolGraph) -> MolGraph:
    """Copy fixture data for an independent render panel."""
    return deepcopy(graph)


def graph_snapshot(graph: MolGraph) -> tuple[tuple[tuple[str, object], ...], ...]:
    """Return complete stable atom/bond state, including coordinates/stereo."""
    atom_rows = tuple(
        tuple(sorted(vars(atom).items()))
        for atom in sorted(graph.atoms.values(), key=lambda value: value.id)
    )
    bond_rows = tuple(
        tuple(sorted(vars(bond).items()))
        for bond in sorted(graph.bonds.values(), key=lambda value: value.id)
    )
    return atom_rows, bond_rows
