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


class BlockKind(str, Enum):
    AROMATIC_RING = "aromatic_ring"
    FUSED_SYSTEM = "fused_system"
    MACROCYCLE = "macrocycle"
    CYCLOPHANE = "cyclophane"
    LINKER = "linker"
    TERMINAL_SUBSTITUENT = "terminal_substituent"
    STEREO_CENTER = "stereo_center"
    INTRAMOLECULAR_BRIDGE = "intramolecular_bridge"
    INTERNAL_CAVITY = "internal_cavity"


class BlockEdgeKind(str, Enum):
    ATTACHMENT = "attachment"
    LINKER = "linker"
    BRIDGE = "bridge"
    CONTAINS = "contains"


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
class BlockNode:
    id: int
    kind: BlockKind
    atom_ids: frozenset[int]
    label: str = ""
    anchor_atom_ids: tuple[int, ...] = ()
    motif_ids: tuple[int, ...] = ()
    centroid: tuple[float, float] | None = None
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BlockEdge:
    id: int
    block_ids: tuple[int, int]
    kind: BlockEdgeKind
    atom_ids: tuple[int, ...] = ()
    weight: float = 1.0
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BlockGraph:
    atom_ids: frozenset[int]
    blocks: tuple[BlockNode, ...] = ()
    edges: tuple[BlockEdge, ...] = ()
    atom_to_block_ids: dict[int, tuple[int, ...]] = field(default_factory=dict)

    @classmethod
    def from_motif_graph(
        cls,
        graph: MolGraph,
        motif_graph: MotifGraph,
        interaction_graph: InteractionGraph | None = None,
    ) -> "BlockGraph":
        selected = motif_graph.atom_ids
        covalent_bonds = _covalent_bonds(graph, selected)
        adjacency = _adjacency(selected, covalent_bonds)
        rings = [motif for motif in motif_graph.motifs if motif.kind == MotifKind.RING]
        blocks: list[BlockNode] = []

        fused_groups = _fused_ring_groups(rings)
        fused_ring_ids = {motif_id for group in fused_groups for motif_id in group if len(group) > 1}
        for group in fused_groups:
            if len(group) <= 1:
                continue
            motifs = [motif for motif in rings if motif.id in group]
            atoms = frozenset().union(*(motif.atom_ids for motif in motifs))
            blocks.append(
                _block(
                    BlockKind.FUSED_SYSTEM,
                    graph,
                    atoms,
                    "fused_system",
                    motif_ids=tuple(sorted(group)),
                    metadata={"ring_count": len(group)},
                )
            )
            if not any(bool(motif.metadata.get("aromatic", False)) for motif in motifs):
                shared_atoms = frozenset.intersection(*(motif.atom_ids for motif in motifs))
                if len(shared_atoms) >= 2:
                    anchors = tuple(
                        sorted(
                            atom_id
                            for atom_id in shared_atoms
                            if len(adjacency.get(atom_id, set()) - shared_atoms) >= 2
                        )
                    )
                    blocks.append(
                        _block(
                            BlockKind.INTRAMOLECULAR_BRIDGE,
                            graph,
                            shared_atoms,
                            "intramolecular_bridge",
                            anchor_atom_ids=anchors,
                            motif_ids=tuple(sorted(group)),
                            metadata={"source": "fused_cycle_intersection"},
                        )
                    )

        for motif in rings:
            size = int(motif.metadata.get("size", len(motif.atom_ids)) or len(motif.atom_ids))
            aromatic = bool(motif.metadata.get("aromatic", False))
            if aromatic:
                blocks.append(
                    _block(
                        BlockKind.AROMATIC_RING,
                        graph,
                        motif.atom_ids,
                        "aromatic_ring",
                        motif_ids=(motif.id,),
                        metadata={"fused": motif.id in fused_ring_ids, "size": size},
                    )
                )
            if size >= 12:
                kind = BlockKind.CYCLOPHANE if _ring_has_aromatic_bridge(motif.atom_ids, rings) else BlockKind.MACROCYCLE
                blocks.append(
                    _block(kind, graph, motif.atom_ids, kind.value, motif_ids=(motif.id,), metadata={"size": size})
                )

        for component in _components(selected, covalent_bonds):
            core = _cycle_core_atoms(component, adjacency)
            aromatic_rings_in_core = sum(
                1
                for ring in rings
                if bool(ring.metadata.get("aromatic", False)) and len(ring.atom_ids & core) >= 3
            )
            already_macro = any(
                block.kind in {BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
                and len(block.atom_ids & core) >= 12
                for block in blocks
            )
            if len(core) >= 12 and aromatic_rings_in_core == 0 and not already_macro:
                blocks.append(
                    _block(
                        BlockKind.MACROCYCLE,
                        graph,
                        core,
                        "macrocycle",
                        metadata={"size": len(core), "source": "cycle_core"},
                    )
                )

        rigid_atoms = frozenset(
            atom_id
            for block in blocks
            if block.kind in {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
            for atom_id in block.atom_ids
        )
        rigid_block_atoms = [
            block.atom_ids
            for block in blocks
            if block.kind in {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
        ]

        for path, endpoints in _linker_paths_between_blocks(selected, adjacency, rigid_block_atoms):
            blocks.append(
                _block(
                    BlockKind.LINKER,
                    graph,
                    path,
                    "linker",
                    anchor_atom_ids=endpoints,
                    metadata={"length": len(path)},
                )
            )

        for atoms, anchors in _terminal_substituents(selected, adjacency, rigid_atoms):
            blocks.append(
                _block(
                    BlockKind.TERMINAL_SUBSTITUENT,
                    graph,
                    atoms,
                    "terminal_substituent",
                    anchor_atom_ids=anchors,
                )
            )

        for atom_id in sorted(selected):
            atom = graph.atoms.get(atom_id)
            if atom is None:
                continue
            stereo_bonds = [
                bond.id
                for bond in covalent_bonds
                if atom_id in {bond.a1_id, bond.a2_id} and bond.stereo != BondStereo.NONE
            ]
            if atom.stereo_cip or stereo_bonds:
                blocks.append(
                    _block(
                        BlockKind.STEREO_CENTER,
                        graph,
                        {atom_id},
                        "stereo_center",
                        anchor_atom_ids=(atom_id,),
                        metadata={"stereo_cip": atom.stereo_cip, "stereo_bond_ids": tuple(stereo_bonds)},
                    )
                )

        for atoms, anchors in _intramolecular_bridges(selected, adjacency, rigid_block_atoms):
            blocks.append(
                _block(
                    BlockKind.INTRAMOLECULAR_BRIDGE,
                    graph,
                    atoms,
                    "intramolecular_bridge",
                    anchor_atom_ids=anchors,
                )
            )

        for motif in motif_graph.motifs:
            if motif.kind == MotifKind.CAVITY:
                blocks.append(_block(BlockKind.INTERNAL_CAVITY, graph, motif.atom_ids, motif.label, motif_ids=(motif.id,)))
        for block in tuple(blocks):
            if block.kind in {BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}:
                blocks.append(
                    _block(
                        BlockKind.INTERNAL_CAVITY,
                        graph,
                        block.atom_ids,
                        "macrocycle_cavity",
                        metadata={"source_block_kind": block.kind.value},
                    )
                )

        fixed_blocks: list[BlockNode] = []
        atom_to_blocks: dict[int, list[int]] = {atom_id: [] for atom_id in selected}
        for idx, block in enumerate(blocks, start=1):
            fixed = BlockNode(
                id=idx,
                kind=block.kind,
                atom_ids=block.atom_ids,
                label=block.label,
                anchor_atom_ids=block.anchor_atom_ids,
                motif_ids=block.motif_ids,
                centroid=block.centroid,
                metadata=block.metadata,
            )
            fixed_blocks.append(fixed)
            for atom_id in fixed.atom_ids:
                atom_to_blocks.setdefault(atom_id, []).append(fixed.id)

        edges = _block_edges(fixed_blocks, covalent_bonds, interaction_graph)
        return cls(
            atom_ids=selected,
            blocks=tuple(fixed_blocks),
            edges=tuple(edges),
            atom_to_block_ids={atom_id: tuple(ids) for atom_id, ids in atom_to_blocks.items()},
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
    block_graph: BlockGraph
    layout_constraint_graph: LayoutConstraintGraph


def build_multilayer_chemical_graph(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
) -> MultilayerChemicalGraph:
    selected = frozenset(_normalize_atom_ids(graph, atom_ids))
    interactions = InteractionGraph.from_molgraph(graph, selected)
    motifs = MotifGraph.from_molgraph(graph, selected, interactions)
    blocks = BlockGraph.from_motif_graph(graph, motifs, interactions)
    constraints = LayoutConstraintGraph.from_layers(graph, selected, interactions, motifs)
    return MultilayerChemicalGraph(
        mol_graph=graph,
        atom_ids=selected,
        interaction_graph=interactions,
        motif_graph=motifs,
        block_graph=blocks,
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
    adjacency = _adjacency(selected, bonds)
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


def _adjacency(atom_ids: Iterable[int], bonds: Iterable[Bond]) -> dict[int, set[int]]:
    selected = set(atom_ids)
    adjacency = {atom_id: set() for atom_id in selected}
    for bond in bonds:
        if bond.a1_id in selected and bond.a2_id in selected:
            adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
            adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    return adjacency


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


def _block(
    kind: BlockKind,
    graph: MolGraph,
    atom_ids: Iterable[int],
    label: str,
    *,
    anchor_atom_ids: tuple[int, ...] = (),
    motif_ids: tuple[int, ...] = (),
    metadata: dict[str, object] | None = None,
) -> BlockNode:
    atoms = frozenset(atom_ids)
    return BlockNode(
        id=0,
        kind=kind,
        atom_ids=atoms,
        label=label,
        anchor_atom_ids=anchor_atom_ids,
        motif_ids=motif_ids,
        centroid=_centroid(graph, atoms),
        metadata=dict(metadata or {}),
    )


def _fused_ring_groups(rings: list[MotifNode]) -> list[set[int]]:
    ring_by_id = {ring.id: ring for ring in rings}
    adjacency = {ring.id: set() for ring in rings}
    for idx, left in enumerate(rings):
        for right in rings[idx + 1 :]:
            if len(left.atom_ids & right.atom_ids) >= 2:
                adjacency[left.id].add(right.id)
                adjacency[right.id].add(left.id)
    groups: list[set[int]] = []
    seen: set[int] = set()
    for ring_id in sorted(ring_by_id):
        if ring_id in seen:
            continue
        stack = [ring_id]
        group: set[int] = set()
        while stack:
            current = stack.pop()
            if current in group:
                continue
            group.add(current)
            stack.extend(sorted(adjacency.get(current, set()) - group))
        seen.update(group)
        groups.append(group)
    return groups


def _ring_has_aromatic_bridge(atom_ids: frozenset[int], rings: list[MotifNode]) -> bool:
    return any(
        bool(ring.metadata.get("aromatic", False)) and len(atom_ids & ring.atom_ids) >= 2
        for ring in rings
    )


def _cycle_core_atoms(atom_ids: Iterable[int], adjacency: dict[int, set[int]]) -> frozenset[int]:
    core = set(atom_ids)
    degrees = {atom_id: len(adjacency.get(atom_id, set()) & core) for atom_id in core}
    queue = [atom_id for atom_id, degree in degrees.items() if degree <= 1]
    while queue:
        atom_id = queue.pop()
        if atom_id not in core:
            continue
        core.remove(atom_id)
        for neighbor in adjacency.get(atom_id, set()):
            if neighbor not in core:
                continue
            degrees[neighbor] = degrees.get(neighbor, 0) - 1
            if degrees[neighbor] <= 1:
                queue.append(neighbor)
    return frozenset(core)


def _linker_paths_between_blocks(
    atom_ids: Iterable[int],
    adjacency: dict[int, set[int]],
    block_atoms: list[frozenset[int]],
) -> list[tuple[frozenset[int], tuple[int, ...]]]:
    selected = set(atom_ids)
    atom_to_blocks: dict[int, set[int]] = {atom_id: set() for atom_id in selected}
    for idx, atoms in enumerate(block_atoms):
        for atom_id in atoms:
            atom_to_blocks.setdefault(atom_id, set()).add(idx)
    block_union = set().union(*block_atoms) if block_atoms else set()
    paths: list[tuple[frozenset[int], tuple[int, ...]]] = []
    seen: set[tuple[int, ...]] = set()
    for start in sorted(block_union):
        for neighbor in sorted(adjacency.get(start, set()) - block_union):
            path = [start, neighbor]
            previous = start
            current = neighbor
            while True:
                next_block_neighbors = sorted(adjacency.get(current, set()) & block_union - {start})
                free_neighbors = sorted((adjacency.get(current, set()) - block_union) - {previous})
                if next_block_neighbors:
                    end = next_block_neighbors[0]
                    if atom_to_blocks.get(start, set()).isdisjoint(atom_to_blocks.get(end, set())):
                        full_path = tuple(path + [end])
                        key = tuple(sorted(full_path))
                        if key not in seen:
                            seen.add(key)
                            paths.append((frozenset(full_path), (start, end)))
                    break
                if len(free_neighbors) != 1:
                    break
                previous, current = current, free_neighbors[0]
                if current in path:
                    break
                path.append(current)
    return paths


def _terminal_substituents(
    atom_ids: Iterable[int],
    adjacency: dict[int, set[int]],
    rigid_atoms: frozenset[int],
) -> list[tuple[frozenset[int], tuple[int, ...]]]:
    free_atoms = set(atom_ids) - set(rigid_atoms)
    out: list[tuple[frozenset[int], tuple[int, ...]]] = []
    seen: set[int] = set()
    for atom_id in sorted(free_atoms):
        if atom_id in seen:
            continue
        stack = [atom_id]
        component: set[int] = set()
        anchors: set[int] = set()
        while stack:
            current = stack.pop()
            if current in component or current not in free_atoms:
                continue
            component.add(current)
            for neighbor in adjacency.get(current, set()):
                if neighbor in rigid_atoms:
                    anchors.add(neighbor)
                elif neighbor not in component:
                    stack.append(neighbor)
        seen.update(component)
        if component and len(anchors) == 1:
            out.append((frozenset(component | anchors), tuple(sorted(anchors))))
    return out


def _intramolecular_bridges(
    atom_ids: Iterable[int],
    adjacency: dict[int, set[int]],
    block_atoms: list[frozenset[int]],
) -> list[tuple[frozenset[int], tuple[int, ...]]]:
    selected = set(atom_ids)
    block_union = set().union(*block_atoms) if block_atoms else set()
    out: list[tuple[frozenset[int], tuple[int, ...]]] = []
    seen: set[tuple[int, ...]] = set()
    for block in block_atoms:
        for start in sorted(block):
            for neighbor in sorted((adjacency.get(start, set()) - block) & selected):
                if neighbor in block_union:
                    continue
                path = [start, neighbor]
                previous = start
                current = neighbor
                while True:
                    block_neighbors = sorted((adjacency.get(current, set()) & block) - {start})
                    if block_neighbors:
                        end = block_neighbors[0]
                        key = tuple(sorted(path + [end]))
                        if key not in seen:
                            seen.add(key)
                            out.append((frozenset(path + [end]), (start, end)))
                        break
                    next_free = sorted((adjacency.get(current, set()) - block_union) - {previous})
                    if len(next_free) != 1:
                        break
                    previous, current = current, next_free[0]
                    if current in path:
                        break
                    path.append(current)
    return out


def _block_edges(
    blocks: list[BlockNode],
    bonds: list[Bond],
    interaction_graph: InteractionGraph | None,
) -> list[BlockEdge]:
    atom_to_blocks: dict[int, set[int]] = {}
    for block in blocks:
        for atom_id in block.atom_ids:
            atom_to_blocks.setdefault(atom_id, set()).add(block.id)
    edges: list[BlockEdge] = []
    seen: set[tuple[int, int, str, tuple[int, ...]]] = set()
    for bond in bonds:
        for left in atom_to_blocks.get(bond.a1_id, set()):
            for right in atom_to_blocks.get(bond.a2_id, set()):
                if left == right:
                    continue
                left_block = blocks[left - 1]
                right_block = blocks[right - 1]
                kind = BlockEdgeKind.LINKER if BlockKind.LINKER in {left_block.kind, right_block.kind} else BlockEdgeKind.ATTACHMENT
                if BlockKind.INTRAMOLECULAR_BRIDGE in {left_block.kind, right_block.kind}:
                    kind = BlockEdgeKind.BRIDGE
                key = (min(left, right), max(left, right), kind.value, tuple(sorted((bond.a1_id, bond.a2_id))))
                if key in seen:
                    continue
                seen.add(key)
                edges.append(
                    BlockEdge(
                        id=len(edges) + 1,
                        block_ids=(min(left, right), max(left, right)),
                        kind=kind,
                        atom_ids=(bond.a1_id, bond.a2_id),
                    )
                )
    if interaction_graph is not None:
        for edge in interaction_graph.edges:
            left_ids = atom_to_blocks.get(edge.atom_ids[0], set())
            right_ids = atom_to_blocks.get(edge.atom_ids[1], set())
            for left in left_ids:
                for right in right_ids:
                    if left == right:
                        continue
                    key = (min(left, right), max(left, right), BlockEdgeKind.BRIDGE.value, edge.atom_ids)
                    if key in seen:
                        continue
                    seen.add(key)
                    edges.append(
                        BlockEdge(
                            id=len(edges) + 1,
                            block_ids=(min(left, right), max(left, right)),
                            kind=BlockEdgeKind.BRIDGE,
                            atom_ids=edge.atom_ids,
                            weight=edge.strength,
                            metadata={"interaction_kind": edge.kind.value},
                        )
                    )
    return edges


def _centroid(graph: MolGraph, atom_ids: Iterable[int]) -> tuple[float, float] | None:
    coords = [(graph.atoms[atom_id].x, graph.atoms[atom_id].y) for atom_id in atom_ids if atom_id in graph.atoms]
    if not coords:
        return None
    return (
        sum(float(x) for x, _ in coords) / len(coords),
        sum(float(y) for _, y in coords) / len(coords),
    )
