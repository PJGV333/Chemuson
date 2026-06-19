from __future__ import annotations

"""Derived multilayer chemical graphs.

The mutable ``MolGraph`` remains the source of truth for atoms and covalent
connectivity.  This module builds read-only semantic layers used by depiction,
validation and future supramolecular tooling without adding non-covalent edges
to the covalent graph.
"""

from collections import deque
from dataclasses import dataclass, field
from enum import Enum
import math
from typing import Iterable

from chemuson.core.model import (
    Bond,
    BondStereo,
    BondStyle,
    MolGraph,
    bond_affects_valence,
    normalize_bond_style,
)


class InteractionKind(str, Enum):
    COORDINATION = "coordination"
    HYDROGEN_BOND = "hydrogen_bond"
    PI_PI = "pi_pi"
    CATION_PI = "cation_pi"
    HOST_GUEST = "host_guest"
    ION_PAIR = "ion_pair"
    NONCOVALENT = "noncovalent"


class MotifKind(str, Enum):
    COVALENT_COMPONENT = "covalent_component"
    RING = "ring"
    FUNCTIONAL_GROUP = "functional_group"
    LIGAND = "ligand"
    CAVITY = "cavity"
    CENTROID = "centroid"
    RIGID_BLOCK = "rigid_block"


class LayoutConstraintKind(str, Enum):
    COVALENT_DISTANCE = "covalent_distance"
    INTERACTION_DISTANCE = "interaction_distance"
    STEREO_ORIENTATION = "stereo_orientation"
    MOTIF_RIGID = "motif_rigid"
    MOTIF_CENTROID = "motif_centroid"
    MOTIF_CONTAINMENT = "motif_containment"


@dataclass(frozen=True)
class InteractionEdge:
    id: int
    atom_ids: tuple[int, int]
    kind: InteractionKind
    source_bond_id: int | None = None
    donor_atom_id: int | None = None
    strength: float = 1.0
    preferred_distance: float | None = None
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class InteractionGraph:
    atom_ids: frozenset[int]
    edges: tuple[InteractionEdge, ...] = ()

    @classmethod
    def from_molgraph(
        cls,
        graph: MolGraph,
        atom_ids: Iterable[int] | None = None,
    ) -> "InteractionGraph":
        selected = frozenset(_normalize_atom_ids(graph, atom_ids))
        edges: list[InteractionEdge] = []
        for bond in sorted(graph.bonds.values(), key=lambda item: item.id):
            if bond.a1_id not in selected or bond.a2_id not in selected:
                continue
            style = normalize_bond_style(bond.style)
            if style not in {BondStyle.COORDINATION, BondStyle.INTERACTION}:
                continue
            kind = _interaction_kind_for_bond(graph, bond)
            edges.append(
                InteractionEdge(
                    id=len(edges) + 1,
                    atom_ids=(bond.a1_id, bond.a2_id),
                    kind=kind,
                    source_bond_id=bond.id,
                    donor_atom_id=bond.donor_atom_id,
                    strength=_interaction_strength(kind),
                    preferred_distance=bond.length_px or _preferred_interaction_distance(kind),
                    metadata={"style": style.value},
                )
            )
        return cls(atom_ids=selected, edges=tuple(edges))


@dataclass(frozen=True)
class MotifNode:
    id: int
    kind: MotifKind
    atom_ids: frozenset[int]
    label: str = ""
    anchor_atom_ids: tuple[int, ...] = ()
    centroid: tuple[float, float] | None = None
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class MotifGraph:
    atom_ids: frozenset[int]
    motifs: tuple[MotifNode, ...] = ()
    atom_to_motif_ids: dict[int, tuple[int, ...]] = field(default_factory=dict)

    @classmethod
    def from_molgraph(
        cls,
        graph: MolGraph,
        atom_ids: Iterable[int] | None = None,
        interaction_graph: InteractionGraph | None = None,
    ) -> "MotifGraph":
        selected = frozenset(_normalize_atom_ids(graph, atom_ids))
        covalent_bonds = _covalent_bonds(graph, selected)
        motifs: list[MotifNode] = []

        for component in _components(selected, covalent_bonds):
            motifs.append(_motif(MotifKind.COVALENT_COMPONENT, graph, component, "component"))

        rings = _cycle_basis(selected, covalent_bonds, max_size=18)
        for ring in rings:
            aromatic = all(
                bond.is_aromatic
                for bond in _ring_bonds(ring, covalent_bonds)
            )
            motifs.append(
                _motif(
                    MotifKind.RING,
                    graph,
                    ring,
                    "aromatic_ring" if aromatic else f"{len(ring)}-ring",
                    metadata={"aromatic": aromatic, "size": len(ring)},
                )
            )
            motifs.append(_motif(MotifKind.CENTROID, graph, ring, "ring_centroid"))
            motifs.append(_motif(MotifKind.RIGID_BLOCK, graph, ring, "ring_rigid_block"))

        for atom_id in sorted(selected):
            atom = graph.atoms.get(atom_id)
            if atom is None:
                continue
            if atom.element not in {"C", "H"} or atom.charge != 0 or atom.group_h_cap is not None:
                neighbors = _neighbors(atom_id, covalent_bonds) & selected
                motif_atoms = frozenset({atom_id, *neighbors})
                if len(motif_atoms) > 1:
                    motifs.append(
                        _motif(
                            MotifKind.FUNCTIONAL_GROUP,
                            graph,
                            motif_atoms,
                            atom.element,
                            anchor_atom_ids=(atom_id,),
                        )
                    )

        interactions = interaction_graph or InteractionGraph.from_molgraph(graph, selected)
        covalent_components = _components(selected, covalent_bonds)
        atom_component = {
            atom_id: component
            for component in covalent_components
            for atom_id in component
        }
        for edge in interactions.edges:
            if edge.kind == InteractionKind.COORDINATION:
                for atom_id in edge.atom_ids:
                    atom = graph.atoms.get(atom_id)
                    if atom and not atom.is_coordination_center:
                        motifs.append(
                            _motif(
                                MotifKind.LIGAND,
                                graph,
                                atom_component.get(atom_id, frozenset({atom_id})),
                                "ligand",
                                anchor_atom_ids=(atom_id,),
                            )
                        )
            if edge.kind == InteractionKind.HOST_GUEST:
                motifs.append(
                    _motif(MotifKind.CAVITY, graph, frozenset(edge.atom_ids), "host_guest_cavity")
                )

        atom_to_motif: dict[int, list[int]] = {atom_id: [] for atom_id in selected}
        for idx, motif in enumerate(motifs, start=1):
            fixed = MotifNode(
                id=idx,
                kind=motif.kind,
                atom_ids=motif.atom_ids,
                label=motif.label,
                anchor_atom_ids=motif.anchor_atom_ids,
                centroid=motif.centroid,
                metadata=motif.metadata,
            )
            motifs[idx - 1] = fixed
            for atom_id in fixed.atom_ids:
                atom_to_motif.setdefault(atom_id, []).append(fixed.id)
        return cls(
            atom_ids=selected,
            motifs=tuple(motifs),
            atom_to_motif_ids={atom_id: tuple(ids) for atom_id, ids in atom_to_motif.items()},
        )


@dataclass(frozen=True)
class LayoutConstraint:
    id: int
    kind: LayoutConstraintKind
    atom_ids: tuple[int, ...] = ()
    motif_ids: tuple[int, ...] = ()
    target_distance: float | None = None
    target_angle_deg: float | None = None
    weight: float = 1.0
    source: str = ""
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class LayoutConstraintGraph:
    atom_ids: frozenset[int]
    constraints: tuple[LayoutConstraint, ...] = ()

    @classmethod
    def from_layers(
        cls,
        graph: MolGraph,
        atom_ids: Iterable[int] | None = None,
        interaction_graph: InteractionGraph | None = None,
        motif_graph: MotifGraph | None = None,
    ) -> "LayoutConstraintGraph":
        selected = frozenset(_normalize_atom_ids(graph, atom_ids))
        interactions = interaction_graph or InteractionGraph.from_molgraph(graph, selected)
        motifs = motif_graph or MotifGraph.from_molgraph(graph, selected, interactions)
        constraints: list[LayoutConstraint] = []

        for bond in _covalent_bonds(graph, selected):
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.COVALENT_DISTANCE,
                    atom_ids=(bond.a1_id, bond.a2_id),
                    target_distance=bond.length_px,
                    weight=1.0,
                    source=f"bond:{bond.id}",
                )
            )
            if bond.stereo != BondStereo.NONE or bond.stereo_ez:
                constraints.append(
                    LayoutConstraint(
                        id=len(constraints) + 1,
                        kind=LayoutConstraintKind.STEREO_ORIENTATION,
                        atom_ids=(bond.a1_id, bond.a2_id),
                        weight=3.0,
                        source=f"stereo:{bond.id}",
                        metadata={"stereo": bond.stereo.value, "stereo_ez": bond.stereo_ez},
                    )
                )

        for edge in interactions.edges:
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.INTERACTION_DISTANCE,
                    atom_ids=edge.atom_ids,
                    target_distance=edge.preferred_distance,
                    weight=edge.strength,
                    source=f"interaction:{edge.id}",
                    metadata={"kind": edge.kind.value, "source_bond_id": edge.source_bond_id},
                )
            )

        for motif in motifs.motifs:
            if motif.kind == MotifKind.RIGID_BLOCK:
                constraints.append(
                    LayoutConstraint(
                        id=len(constraints) + 1,
                        kind=LayoutConstraintKind.MOTIF_RIGID,
                        atom_ids=tuple(sorted(motif.atom_ids)),
                        motif_ids=(motif.id,),
                        weight=4.0,
                        source=f"motif:{motif.id}",
                    )
                )
            elif motif.kind == MotifKind.CENTROID:
                constraints.append(
                    LayoutConstraint(
                        id=len(constraints) + 1,
                        kind=LayoutConstraintKind.MOTIF_CENTROID,
                        atom_ids=tuple(sorted(motif.atom_ids)),
                        motif_ids=(motif.id,),
                        weight=2.0,
                        source=f"motif:{motif.id}",
                    )
                )
            elif motif.kind == MotifKind.CAVITY:
                constraints.append(
                    LayoutConstraint(
                        id=len(constraints) + 1,
                        kind=LayoutConstraintKind.MOTIF_CONTAINMENT,
                        atom_ids=tuple(sorted(motif.atom_ids)),
                        motif_ids=(motif.id,),
                        weight=1.5,
                        source=f"motif:{motif.id}",
                    )
                )
        return cls(atom_ids=selected, constraints=tuple(constraints))


@dataclass(frozen=True)
class MultilayerChemicalGraph:
    mol_graph: MolGraph
    atom_ids: frozenset[int]
    interaction_graph: InteractionGraph
    motif_graph: MotifGraph
    layout_constraint_graph: LayoutConstraintGraph


def build_multilayer_chemical_graph(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
) -> MultilayerChemicalGraph:
    selected = frozenset(_normalize_atom_ids(graph, atom_ids))
    interactions = InteractionGraph.from_molgraph(graph, selected)
    motifs = MotifGraph.from_molgraph(graph, selected, interactions)
    constraints = LayoutConstraintGraph.from_layers(graph, selected, interactions, motifs)
    return MultilayerChemicalGraph(
        mol_graph=graph,
        atom_ids=selected,
        interaction_graph=interactions,
        motif_graph=motifs,
        layout_constraint_graph=constraints,
    )


def _normalize_atom_ids(graph: MolGraph, atom_ids: Iterable[int] | None) -> set[int]:
    if atom_ids is None:
        return set(graph.atoms.keys())
    return {int(atom_id) for atom_id in atom_ids if int(atom_id) in graph.atoms}


def _covalent_bonds(graph: MolGraph, atom_ids: Iterable[int]) -> list[Bond]:
    selected = set(atom_ids)
    return [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    ]


def _interaction_kind_for_bond(graph: MolGraph, bond: Bond) -> InteractionKind:
    style = normalize_bond_style(bond.style)
    explicit_kind = getattr(bond, "interaction_kind", None)
    if explicit_kind:
        try:
            return InteractionKind(str(explicit_kind).strip().lower())
        except ValueError:
            pass
    if style == BondStyle.COORDINATION:
        return InteractionKind.COORDINATION
    atoms = [graph.atoms.get(bond.a1_id), graph.atoms.get(bond.a2_id)]
    elements = {atom.element for atom in atoms if atom is not None}
    charges = [atom.charge for atom in atoms if atom is not None]
    if {"H"} & elements and elements & {"N", "O", "F", "S", "Cl"}:
        return InteractionKind.HYDROGEN_BOND
    if any(charge > 0 for charge in charges) and any(atom and atom.element == "C" for atom in atoms):
        return InteractionKind.CATION_PI
    if any(charge > 0 for charge in charges) and any(charge < 0 for charge in charges):
        return InteractionKind.ION_PAIR
    return InteractionKind.NONCOVALENT


def _interaction_strength(kind: InteractionKind) -> float:
    return {
        InteractionKind.COORDINATION: 2.8,
        InteractionKind.ION_PAIR: 2.4,
        InteractionKind.HYDROGEN_BOND: 1.8,
        InteractionKind.CATION_PI: 1.4,
        InteractionKind.PI_PI: 1.2,
        InteractionKind.HOST_GUEST: 1.0,
        InteractionKind.NONCOVALENT: 0.9,
    }[kind]


def _preferred_interaction_distance(kind: InteractionKind) -> float:
    return {
        InteractionKind.COORDINATION: 44.0,
        InteractionKind.ION_PAIR: 58.0,
        InteractionKind.HYDROGEN_BOND: 52.0,
        InteractionKind.CATION_PI: 62.0,
        InteractionKind.PI_PI: 68.0,
        InteractionKind.HOST_GUEST: 76.0,
        InteractionKind.NONCOVALENT: 64.0,
    }[kind]


def _components(atom_ids: Iterable[int], bonds: Iterable[Bond]) -> list[frozenset[int]]:
    selected = set(atom_ids)
    adjacency = {atom_id: set() for atom_id in selected}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
        adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    components: list[frozenset[int]] = []
    seen: set[int] = set()
    for atom_id in sorted(selected):
        if atom_id in seen:
            continue
        queue = deque([atom_id])
        component: set[int] = set()
        while queue:
            current = queue.popleft()
            if current in component:
                continue
            component.add(current)
            queue.extend(sorted(adjacency.get(current, set()) - component))
        seen.update(component)
        components.append(frozenset(component))
    return components


def _neighbors(atom_id: int, bonds: Iterable[Bond]) -> set[int]:
    out: set[int] = set()
    for bond in bonds:
        if bond.a1_id == atom_id:
            out.add(bond.a2_id)
        elif bond.a2_id == atom_id:
            out.add(bond.a1_id)
    return out


def _cycle_basis(atom_ids: Iterable[int], bonds: list[Bond], *, max_size: int) -> list[frozenset[int]]:
    selected = set(atom_ids)
    adjacency = {atom_id: set() for atom_id in selected}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
        adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    rings: list[frozenset[int]] = []
    seen: set[tuple[int, ...]] = set()
    for bond in sorted(bonds, key=lambda item: item.id):
        path = _shortest_path_without_edge(adjacency, bond.a1_id, bond.a2_id)
        if path is None or not (3 <= len(path) <= max_size):
            continue
        key = tuple(sorted(path))
        if key in seen:
            continue
        seen.add(key)
        rings.append(frozenset(path))
    return rings


def _shortest_path_without_edge(adjacency: dict[int, set[int]], start: int, end: int) -> list[int] | None:
    blocked = {start, end}
    queue: deque[tuple[int, list[int]]] = deque([(start, [start])])
    visited: set[int] = set()
    while queue:
        atom_id, path = queue.popleft()
        if atom_id == end and len(path) >= 3:
            return path
        if atom_id in visited:
            continue
        visited.add(atom_id)
        for neighbor in sorted(adjacency.get(atom_id, set())):
            if {atom_id, neighbor} == blocked or neighbor in path:
                continue
            queue.append((neighbor, path + [neighbor]))
    return None


def _ring_bonds(ring: Iterable[int], bonds: Iterable[Bond]) -> list[Bond]:
    atoms = set(ring)
    return [bond for bond in bonds if bond.a1_id in atoms and bond.a2_id in atoms]


def _motif(
    kind: MotifKind,
    graph: MolGraph,
    atom_ids: Iterable[int],
    label: str,
    *,
    anchor_atom_ids: tuple[int, ...] = (),
    metadata: dict[str, object] | None = None,
) -> MotifNode:
    atoms = frozenset(atom_ids)
    centroid = _centroid(graph, atoms)
    return MotifNode(
        id=0,
        kind=kind,
        atom_ids=atoms,
        label=label,
        anchor_atom_ids=anchor_atom_ids,
        centroid=centroid,
        metadata=dict(metadata or {}),
    )


def _centroid(graph: MolGraph, atom_ids: Iterable[int]) -> tuple[float, float] | None:
    coords = [(graph.atoms[atom_id].x, graph.atoms[atom_id].y) for atom_id in atom_ids if atom_id in graph.atoms]
    if not coords:
        return None
    return (
        sum(float(x) for x, _ in coords) / len(coords),
        sum(float(y) for _, y in coords) / len(coords),
    )
