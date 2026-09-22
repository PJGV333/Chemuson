from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

from chemuson.core.layers import (
    BlockKind,
    MultilayerChemicalGraph,
    build_multilayer_chemical_graph,
)
from chemuson.core.model import MolGraph, bond_affects_valence

if TYPE_CHECKING:
    from chemuson.clean2d.engine import Clean2DLayoutQualityReport


@dataclass(frozen=True)
class Clean2DComplexityProfile:
    atom_count: int
    bond_count: int
    ring_count: int
    aromatic_ring_count: int
    fused_system_count: int
    macrocycle_count: int
    cyclophane_count: int
    intramolecular_bridge_count: int
    internal_cavity_count: int
    linker_count: int
    terminal_substituent_count: int
    stereo_center_count: int
    block_count: int
    has_hierarchical_blocks: bool
    has_block_layout_problem: bool
    preserve_only: bool
    preserve_only_required: bool
    local_repair_allowed: bool
    global_redraw_allowed: bool
    internal_route_allowed: bool
    reason: str
    block_counts: dict[str, int]

    @property
    def policy_evidence(self) -> dict[str, object]:
        return {
            "preserve_only_required": self.preserve_only_required,
            "local_repair_allowed": self.local_repair_allowed,
            "global_redraw_allowed": self.global_redraw_allowed,
            "internal_route_allowed": self.internal_route_allowed,
            "reason": self.reason,
            "has_hierarchical_blocks": self.has_hierarchical_blocks,
            "has_block_layout_problem": self.has_block_layout_problem,
            "block_counts": dict(self.block_counts),
        }


def describe_clean2d_topology(layer_model: MultilayerChemicalGraph) -> dict[str, object]:
    """Return stable, JSON-safe evidence from an existing multilayer model."""

    selected = set(layer_model.atom_ids)
    graph = layer_model.mol_graph
    bonds = [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    ]
    components = _connected_components(
        selected,
        (
            (bond.a1_id, bond.a2_id)
            for bond in graph.bonds.values()
            if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
        ),
    )
    blocks = [
        {
            "id": block.id,
            "kind": block.kind.value,
            "atom_ids": sorted(block.atom_ids),
            "anchor_atom_ids": list(block.anchor_atom_ids),
            "motif_ids": list(block.motif_ids),
            "metadata": _stable_json_value(block.metadata),
        }
        for block in sorted(layer_model.block_graph.blocks, key=lambda item: item.id)
    ]
    connectors = [
        {
            "id": edge.id,
            "kind": edge.kind.value,
            "block_ids": list(edge.block_ids),
            "atom_ids": list(edge.atom_ids),
            "weight": _stable_json_value(edge.weight),
            "metadata": _stable_json_value(edge.metadata),
        }
        for edge in sorted(layer_model.block_graph.edges, key=lambda item: item.id)
    ]
    block_kind_counts: dict[str, int] = {}
    for block in blocks:
        kind = str(block["kind"])
        block_kind_counts[kind] = block_kind_counts.get(kind, 0) + 1
    block_kind_counts = dict(sorted(block_kind_counts.items()))
    rigid_block_kinds = {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM}
    semi_rigid_block_kinds = {BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE, BlockKind.INTRAMOLECULAR_BRIDGE}
    rigid_block_ids = sorted(block["id"] for block in blocks if block["kind"] in {kind.value for kind in rigid_block_kinds})
    semi_rigid_block_ids = sorted(
        block["id"] for block in blocks if block["kind"] in {kind.value for kind in semi_rigid_block_kinds}
    )
    ring_motifs = [motif for motif in layer_model.motif_graph.motifs if motif.kind.value == "ring"]
    ring_systems = _ring_systems(ring_motifs, layer_model.block_graph.blocks, bonds)
    block_adjacency = _block_adjacency(connectors)
    articulation_atom_ids, articulation_bond_ids = _articulation_structure(selected, bonds)
    degree = _bond_adjacency(selected, bonds)
    attachment_atom_ids = sorted({atom_id for edge in connectors for atom_id in edge["atom_ids"]})
    flexible_connectors = [
        {
            "id": edge["id"],
            "block_ids": list(edge["block_ids"]),
            "attachment_atom_ids": list(edge["atom_ids"]),
            "rotatable_bond_ids": _rotatable_bond_ids(graph, edge, layer_model.block_graph.blocks),
        }
        for edge in connectors
        if edge["kind"] == "linker"
    ]
    fused_systems = [system for system in ring_systems if system["kind"] == "fused"]
    spiro_systems = [system for system in ring_systems if system["kind"] == "spiro"]
    bridged_systems = [system for system in ring_systems if system["kind"] == "bridged"]
    macrocycle_blocks = sorted(
        block["id"] for block in blocks if block["kind"] == BlockKind.MACROCYCLE.value
    )
    ring_count = sum(
        1
        for motif in layer_model.motif_graph.motifs
        if motif.kind.value == "ring" and motif.label != "ring_centroid"
    )
    bond_count = sum(
        1
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    )
    return {
        "atom_ids": sorted(selected),
        "atom_count": len(selected),
        "bond_count": bond_count,
        "ring_count": ring_count,
        "component_count": len(components),
        "components": [list(component) for component in components],
        "block_count": len(blocks),
        "block_kind_counts": block_kind_counts,
        "blocks": blocks,
        "connector_count": len(connectors),
        "connectors": connectors,
        "ring_systems": ring_systems,
        "rigid_block_ids": rigid_block_ids,
        "semi_rigid_block_ids": semi_rigid_block_ids,
        "block_adjacency": block_adjacency,
        "flexible_connectors": flexible_connectors,
        "attachment_atom_ids": attachment_atom_ids,
        "articulation_atom_ids": articulation_atom_ids,
        "articulation_bond_ids": articulation_bond_ids,
        "branch_point_atom_ids": sorted(atom_id for atom_id, neighbors in degree.items() if len(neighbors) > 2),
        "fused_systems": fused_systems,
        "spiro_systems": spiro_systems,
        "bridged_systems": bridged_systems,
        "macrocycle_blocks": macrocycle_blocks,
    }


def describe_clean2d_rigid_systems(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
) -> dict[str, object]:
    """Describe rigid ring systems and local attachments deterministically.

    The descriptor is deliberately observational: it derives every member and
    relationship from the existing multilayer topology and never identifies a
    molecule by name or fixture.
    """
    selected = _normalize_atom_ids(graph, atom_ids)
    layer_model = build_multilayer_chemical_graph(graph, selected)
    topology = describe_clean2d_topology(layer_model)
    selected_bonds = [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    ]
    adjacency = _bond_adjacency(selected, selected_bonds)
    motifs = {
        motif.id: motif
        for motif in layer_model.motif_graph.motifs
        if motif.kind.value == "ring"
    }
    bond_by_pair = {
        frozenset((bond.a1_id, bond.a2_id)): bond
        for bond in selected_bonds
    }
    systems: list[dict[str, object]] = []
    for system in cast(list[dict[str, object]], topology["ring_systems"]):
        ring_ids = sorted(cast(list[int], system["ring_ids"]))
        ring_records = [motifs[ring_id] for ring_id in ring_ids if ring_id in motifs]
        system_atoms = set(cast(list[int], system["atom_ids"]))
        shared_atoms = set(cast(list[int], system["shared_atom_ids"]))
        ring_bond_ids: set[int] = set()
        aromatic_flags: list[bool] = []
        for ring in ring_records:
            ring_atoms = set(ring.atom_ids)
            aromatic_flags.append(bool(ring.metadata.get("aromatic", False)))
            for pair, bond in bond_by_pair.items():
                if pair <= ring_atoms:
                    ring_bond_ids.add(bond.id)
        shared_bond_ids = sorted(
            bond.id
            for bond in selected_bonds
            if bond.a1_id in shared_atoms and bond.a2_id in shared_atoms
        )
        attachment_atom_ids = sorted(
            atom_id
            for atom_id in system_atoms
            if any(neighbor not in system_atoms for neighbor in adjacency.get(atom_id, set()))
        )
        external_neighbors = sorted(
            neighbor
            for atom_id in attachment_atom_ids
            for neighbor in adjacency.get(atom_id, set())
            if neighbor not in system_atoms
        )
        centroid = _centroid_for_atoms(graph, system_atoms)
        principal_orientation = _principal_orientation(graph, system_atoms, centroid)
        attachment_vectors = []
        for atom_id in attachment_atom_ids:
            ax, ay = graph.atoms[atom_id].x, graph.atoms[atom_id].y
            for neighbor in sorted(adjacency.get(atom_id, set())):
                if neighbor in system_atoms:
                    continue
                nx, ny = graph.atoms[neighbor].x, graph.atoms[neighbor].y
                dx, dy = nx - ax, ny - ay
                length = math.hypot(dx, dy)
                if length <= 1e-9:
                    vector = (0.0, 0.0)
                    angle = 0.0
                else:
                    vector = (dx / length, dy / length)
                    angle = math.atan2(dy, dx)
                attachment_vectors.append(
                    {
                        "atom_id": atom_id,
                        "neighbor_id": neighbor,
                        "vector": vector,
                        "angle_deg": math.degrees(angle),
                        "length": length,
                    }
                )
        attachment_vectors.sort(key=lambda item: (item["atom_id"], item["neighbor_id"]))
        angles = sorted(float(item["angle_deg"]) for item in attachment_vectors)
        local_congestion = _rigid_attachment_congestion(angles)
        relationships = []
        for left_index, left in enumerate(ring_records):
            for right in ring_records[left_index + 1 :]:
                overlap = sorted(set(left.atom_ids) & set(right.atom_ids))
                shared_pair = frozenset(overlap)
                if len(overlap) == 1:
                    relation = "spiro"
                elif len(overlap) >= 2 and shared_pair in bond_by_pair:
                    relation = "fused"
                elif len(overlap) >= 2:
                    relation = "bridged"
                else:
                    relation = "disjoint"
                relationships.append(
                    {
                        "ring_ids": [left.id, right.id],
                        "shared_atom_ids": overlap,
                        "shared_bond_ids": [bond_by_pair[shared_pair].id] if shared_pair in bond_by_pair else [],
                        "relationship": relation,
                    }
                )
        relationships.sort(key=lambda item: tuple(item["ring_ids"]))
        family = str(system["kind"])
        if len(ring_ids) >= 3 and family in {"fused", "ring_system"}:
            family = "polycyclic"
        systems.append(
            {
                "id": int(cast(int, system["id"])),
                "family": family,
                "kind": str(system["kind"]),
                "atom_ids": sorted(system_atoms),
                "bond_ids": sorted(ring_bond_ids),
                "ring_ids": ring_ids,
                "ring_count": len(ring_ids),
                "shared_atom_ids": sorted(shared_atoms),
                "shared_bond_ids": shared_bond_ids,
                "bridgehead_atom_ids": list(cast(list[int], system["bridgehead_atom_ids"])),
                "attachment_atom_ids": attachment_atom_ids,
                "external_neighbor_ids": external_neighbors,
                "external_substituent_count": len(attachment_vectors),
                "centroid": centroid,
                "principal_orientation": principal_orientation,
                "attachment_vectors": attachment_vectors,
                "local_congestion": local_congestion,
                "aromatic": bool(aromatic_flags) and all(aromatic_flags),
                "ring_aromatic_flags": aromatic_flags,
                "ring_relationships": relationships,
            }
        )
    systems.sort(key=lambda item: int(item["id"]))
    return {
        "version": 1,
        "atom_ids": sorted(selected),
        "multiple_rigid_blocks": len(systems) >= 2,
        "rigid_system_count": len(systems),
        "systems": _stable_json_value(systems),
    }


def _centroid_for_atoms(graph: MolGraph, atom_ids: set[int]) -> tuple[float, float]:
    if not atom_ids:
        return (0.0, 0.0)
    return (
        sum(graph.atoms[atom_id].x for atom_id in atom_ids) / len(atom_ids),
        sum(graph.atoms[atom_id].y for atom_id in atom_ids) / len(atom_ids),
    )


def _principal_orientation(
    graph: MolGraph,
    atom_ids: set[int],
    centroid: tuple[float, float],
) -> tuple[float, float]:
    if len(atom_ids) < 2:
        return (1.0, 0.0)
    xx = yy = xy = 0.0
    for atom_id in sorted(atom_ids):
        dx = graph.atoms[atom_id].x - centroid[0]
        dy = graph.atoms[atom_id].y - centroid[1]
        xx += dx * dx
        yy += dy * dy
        xy += dx * dy
    angle = 0.5 * math.atan2(2.0 * xy, xx - yy) if xx or yy else 0.0
    vector = (math.cos(angle), math.sin(angle))
    if vector[0] < -1e-12 or (abs(vector[0]) <= 1e-12 and vector[1] < 0.0):
        vector = (-vector[0], -vector[1])
    return vector


def _rigid_attachment_congestion(angles: list[float]) -> dict[str, object]:
    if len(angles) < 2:
        return {"attachment_count": len(angles), "minimum_angle_deg": None, "crowded_pair_count": 0}
    gaps = []
    for index, angle in enumerate(angles):
        next_angle = angles[(index + 1) % len(angles)]
        gap = (next_angle - angle) % 360.0
        gaps.append(gap)
    minimum = min(gaps)
    return {
        "attachment_count": len(angles),
        "minimum_angle_deg": minimum,
        "crowded_pair_count": sum(1 for gap in gaps if gap < 45.0),
    }


def plan_clean2d_block_assembly(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
) -> dict[str, object]:
    """Return a deterministic block traversal plan for medium assembly."""
    topology = describe_clean2d_topology(build_multilayer_chemical_graph(graph, atom_ids))
    block_records = cast(list[dict[str, object]], topology["blocks"])
    blocks = {int(cast(int, block["id"])): block for block in block_records}
    rigid_ids = {int(block_id) for block_id in cast(list[int], topology["rigid_block_ids"])}
    semi_rigid_ids = {int(block_id) for block_id in cast(list[int], topology["semi_rigid_block_ids"])}
    structural_ids = rigid_ids | semi_rigid_ids
    requires_global_assembly = len(structural_ids) >= 2 and bool(topology["block_adjacency"])
    if not requires_global_assembly:
        return {
            "version": 1,
            "requires_global_assembly": False,
            "anchor_block_id": None,
            "ordered_block_ids": [],
            "ordered_connector_ids": [],
            "flexible_connector_ids": [],
            "flexible_connector_count": 0,
        }

    anchor_block_id = min(
        structural_ids,
        key=lambda block_id: (-len(cast(list[int], blocks[block_id]["atom_ids"])), block_id),
    )
    connector_records = {
        int(cast(int, connector["id"])): connector
        for connector in cast(list[dict[str, object]], topology["connectors"])
    }
    adjacency: dict[int, list[tuple[int, int, str]]] = {block_id: [] for block_id in blocks}
    block_adjacency = cast(list[dict[str, object]], topology["block_adjacency"])
    for relation in block_adjacency:
        relation_block_ids = tuple(int(block_id) for block_id in cast(list[int], relation["block_ids"]))
        for connector_id in cast(list[int], relation["connector_ids"]):
            connector = connector_records[int(connector_id)]
            kind = str(connector["kind"])
            left, right = relation_block_ids
            adjacency[left].append((right, int(connector_id), kind))
            adjacency[right].append((left, int(connector_id), kind))

    priority = {"attachment": 0, "bridge": 1, "linker": 2, "contains": 3}
    visited = {anchor_block_id}
    queue = [anchor_block_id]
    ordered_block_ids: list[int] = []
    ordered_connector_ids: list[int] = []
    while queue:
        block_id = queue.pop(0)
        ordered_block_ids.append(block_id)
        for neighbor, connector_id, kind in sorted(
            adjacency.get(block_id, ()),
            key=lambda item: (priority.get(item[2], 4), item[1], item[0]),
        ):
            if connector_id not in ordered_connector_ids:
                ordered_connector_ids.append(connector_id)
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    for block_id in sorted(blocks):
        if block_id not in visited:
            ordered_block_ids.append(block_id)
    for connector_id in sorted(connector_records):
        if connector_id not in ordered_connector_ids:
            ordered_connector_ids.append(connector_id)
    flexible_connector_ids = sorted(
        int(cast(int, connector["id"]))
        for connector in cast(list[dict[str, object]], topology["flexible_connectors"])
    )
    return {
        "version": 1,
        "requires_global_assembly": True,
        "anchor_block_id": anchor_block_id,
        "ordered_block_ids": ordered_block_ids,
        "ordered_connector_ids": ordered_connector_ids,
        "flexible_connector_ids": flexible_connector_ids,
        "flexible_connector_count": len(flexible_connector_ids),
    }


def _connected_components(atom_ids: set[int], bonds: Iterable[tuple[int, int]]) -> list[tuple[int, ...]]:
    adjacency = {atom_id: set() for atom_id in atom_ids}
    for left, right in bonds:
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)
    components: list[tuple[int, ...]] = []
    unseen = set(atom_ids)
    while unseen:
        start = min(unseen)
        stack = [start]
        component: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            unseen.discard(atom_id)
            stack.extend(sorted(adjacency.get(atom_id, set()) - component, reverse=True))
        components.append(tuple(sorted(component)))
    return sorted(components)


def _stable_json_value(value: Any) -> object:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _stable_json_value(value[key]) for key in sorted(value, key=str)}
    if isinstance(value, (list, tuple)):
        return [_stable_json_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        normalized = [_stable_json_value(item) for item in value]
        return sorted(normalized, key=lambda item: (type(item).__name__, repr(item)))
    if hasattr(value, "value") and isinstance(value.value, (bool, int, float, str)):
        return value.value
    return type(value).__name__


def _bond_adjacency(atom_ids: set[int], bonds: Iterable[Any]) -> dict[int, set[int]]:
    adjacency = {atom_id: set() for atom_id in atom_ids}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
        adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    return adjacency


def _articulation_structure(atom_ids: set[int], bonds: Iterable[Any]) -> tuple[list[int], list[int]]:
    adjacency: dict[int, list[tuple[int, int]]] = {atom_id: [] for atom_id in atom_ids}
    for bond in bonds:
        adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond.id))
        adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond.id))
    discovery: dict[int, int] = {}
    low: dict[int, int] = {}
    articulation_atoms: set[int] = set()
    articulation_bonds: set[int] = set()
    clock = 0

    def visit(atom_id: int, parent_bond_id: int | None) -> None:
        nonlocal clock
        clock += 1
        discovery[atom_id] = low[atom_id] = clock
        child_count = 0
        for neighbor, bond_id in sorted(adjacency.get(atom_id, []), key=lambda item: (item[0], item[1])):
            if bond_id == parent_bond_id:
                continue
            if neighbor not in discovery:
                child_count += 1
                visit(neighbor, bond_id)
                low[atom_id] = min(low[atom_id], low[neighbor])
                if parent_bond_id is not None and low[neighbor] >= discovery[atom_id]:
                    articulation_atoms.add(atom_id)
                if low[neighbor] > discovery[atom_id]:
                    articulation_bonds.add(bond_id)
            else:
                low[atom_id] = min(low[atom_id], discovery[neighbor])
        if parent_bond_id is None and child_count > 1:
            articulation_atoms.add(atom_id)

    for atom_id in sorted(atom_ids):
        if atom_id not in discovery:
            visit(atom_id, None)
    return sorted(articulation_atoms), sorted(articulation_bonds)


def _ring_systems(ring_motifs: list[Any], blocks: Iterable[Any], bonds: Iterable[Any]) -> list[dict[str, object]]:
    ring_by_id = {motif.id: motif for motif in ring_motifs}
    ring_ids = sorted(ring_by_id)
    adjacency = {ring_id: set() for ring_id in ring_ids}
    for index, left_id in enumerate(ring_ids):
        for right_id in ring_ids[index + 1 :]:
            if ring_by_id[left_id].atom_ids & ring_by_id[right_id].atom_ids:
                adjacency[left_id].add(right_id)
                adjacency[right_id].add(left_id)
    components: list[tuple[int, ...]] = []
    unseen = set(ring_ids)
    while unseen:
        start = min(unseen)
        stack = [start]
        component: set[int] = set()
        while stack:
            ring_id = stack.pop()
            if ring_id in component:
                continue
            component.add(ring_id)
            unseen.discard(ring_id)
            stack.extend(sorted(adjacency[ring_id] - component, reverse=True))
        components.append(tuple(sorted(component)))
    bond_pairs = {frozenset((bond.a1_id, bond.a2_id)) for bond in bonds}
    block_list = list(blocks)
    systems: list[dict[str, object]] = []
    for system_id, component_ids in enumerate(sorted(components), start=1):
        motifs = [ring_by_id[ring_id] for ring_id in component_ids]
        atom_ids = frozenset().union(*(motif.atom_ids for motif in motifs))
        shared_atoms: set[int] = set()
        spiro_atoms: set[int] = set()
        has_fused_overlap = False
        has_bridge_overlap = False
        bridgeheads: set[int] = set()
        for index, left in enumerate(motifs):
            for right in motifs[index + 1 :]:
                overlap = left.atom_ids & right.atom_ids
                if len(overlap) == 1:
                    spiro_atoms.update(overlap)
                if len(overlap) >= 2:
                    shared_atoms.update(overlap)
                    if len(overlap) == 2 and any(frozenset(pair) in bond_pairs for pair in _pairs(sorted(overlap))):
                        has_fused_overlap = True
                    else:
                        has_bridge_overlap = True
                        bridgeheads.update(overlap)
        system_bond_adjacency = _bond_adjacency(set(atom_ids), bonds)
        bridgeheads.update(atom_id for atom_id, neighbors in system_bond_adjacency.items() if len(neighbors) >= 3)
        bridge_blocks = [
            block
            for block in block_list
            if block.kind == BlockKind.INTRAMOLECULAR_BRIDGE and block.atom_ids & atom_ids
        ]
        if bridge_blocks:
            bridgeheads.update(atom_id for block in bridge_blocks for atom_id in block.atom_ids & atom_ids)
        block_ids = sorted(block.id for block in block_list if block.atom_ids & atom_ids)
        if has_bridge_overlap or (bridge_blocks and not has_fused_overlap):
            kind = "bridged"
        elif has_fused_overlap:
            kind = "fused"
        elif spiro_atoms:
            kind = "spiro"
        elif len(component_ids) == 1:
            kind = "monocycle"
        else:
            kind = "ring_system"
        systems.append(
            {
                "id": system_id,
                "kind": kind,
                "ring_ids": list(component_ids),
                "atom_ids": sorted(atom_ids),
                "shared_atom_ids": sorted(shared_atoms | spiro_atoms),
                "bridgehead_atom_ids": sorted(bridgeheads) if kind == "bridged" else [],
                "block_ids": block_ids,
            }
        )
    return systems


def _pairs(values: list[int]) -> Iterable[tuple[int, int]]:
    for index, left in enumerate(values):
        for right in values[index + 1 :]:
            yield left, right


def _block_adjacency(connectors: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, int], dict[str, Any]] = {}
    for connector in connectors:
        block_id_values = cast(list[int], connector["block_ids"])
        if len(block_id_values) != 2:
            continue
        block_ids = (block_id_values[0], block_id_values[1])
        item = grouped.setdefault(block_ids, {"block_ids": list(block_ids), "connector_ids": [], "kinds": []})
        item["connector_ids"].append(cast(int, connector["id"]))
        item["kinds"].append(cast(str, connector["kind"]))
    return [
        {**item, "connector_ids": sorted(item["connector_ids"]), "kinds": sorted(set(item["kinds"]))}
        for _, item in sorted(grouped.items())
    ]


def _rotatable_bond_ids(graph: MolGraph, connector: dict[str, object], blocks: Iterable[Any]) -> list[int]:
    block_ids = set(connector["block_ids"])
    linker_atoms = set().union(*(block.atom_ids for block in blocks if block.id in block_ids and block.kind == BlockKind.LINKER))
    if not linker_atoms:
        return []
    return sorted(
        bond.id
        for bond in graph.bonds.values()
        if bond.a1_id in linker_atoms
        and bond.a2_id in linker_atoms
        and bond.order == 1
        and not bond.is_aromatic
        and bond.stereo.value == "none"
    )


def classify_clean2d_complexity(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    target_bond_length: float = 42.0,
    quality_report: Clean2DLayoutQualityReport | None = None,
    layer_model: MultilayerChemicalGraph | None = None,
) -> Clean2DComplexityProfile:
    del target_bond_length
    selected = _normalize_atom_ids(graph, atom_ids)
    layer_model = layer_model or build_multilayer_chemical_graph(graph, selected)
    block_graph = layer_model.block_graph

    counts: dict[BlockKind, int] = {}
    for block in block_graph.blocks:
        counts[block.kind] = counts.get(block.kind, 0) + 1
    block_counts = {kind.value: count for kind, count in sorted(counts.items(), key=lambda item: item[0].value)}

    ring_count = sum(1 for motif in layer_model.motif_graph.motifs if motif.label != "ring_centroid" and motif.kind.value == "ring")
    aromatic_ring_count = counts.get(BlockKind.AROMATIC_RING, 0)
    fused_system_count = counts.get(BlockKind.FUSED_SYSTEM, 0)
    macrocycle_count = counts.get(BlockKind.MACROCYCLE, 0)
    cyclophane_count = counts.get(BlockKind.CYCLOPHANE, 0)
    intramolecular_bridge_count = counts.get(BlockKind.INTRAMOLECULAR_BRIDGE, 0)
    internal_cavity_count = counts.get(BlockKind.INTERNAL_CAVITY, 0)
    linker_count = counts.get(BlockKind.LINKER, 0) or _fallback_linker_count(graph, selected, block_graph)
    if linker_count:
        block_counts[BlockKind.LINKER.value] = linker_count
    terminal_substituent_count = counts.get(BlockKind.TERMINAL_SUBSTITUENT, 0)
    stereo_center_count = counts.get(BlockKind.STEREO_CENTER, 0)
    atom_count = len(selected)
    bond_count = sum(
        1
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond)
    )

    has_hierarchical_blocks = _has_hierarchical_block_layout_signals(block_graph)
    has_block_layout_problem = _has_intramolecular_block_layout_problem(block_graph, quality_report)
    preserve_only, reason = _preserve_reason(
        atom_count=atom_count,
        ring_count=ring_count,
        aromatic_ring_count=aromatic_ring_count,
        fused_system_count=fused_system_count,
        macrocycle_count=macrocycle_count,
        cyclophane_count=cyclophane_count,
        intramolecular_bridge_count=intramolecular_bridge_count,
        internal_cavity_count=internal_cavity_count,
        stereo_center_count=stereo_center_count,
        has_hierarchical_blocks=has_hierarchical_blocks,
        has_block_layout_problem=has_block_layout_problem,
    )
    high_risk_for_redraw = bool(
        preserve_only
        or has_hierarchical_blocks
        or macrocycle_count
        or cyclophane_count
        or intramolecular_bridge_count
        or internal_cavity_count
    )
    global_redraw_allowed = not high_risk_for_redraw
    local_repair_allowed = False
    internal_route_allowed = high_risk_for_redraw

    return Clean2DComplexityProfile(
        atom_count=atom_count,
        bond_count=bond_count,
        ring_count=ring_count,
        aromatic_ring_count=aromatic_ring_count,
        fused_system_count=fused_system_count,
        macrocycle_count=macrocycle_count,
        cyclophane_count=cyclophane_count,
        intramolecular_bridge_count=intramolecular_bridge_count,
        internal_cavity_count=internal_cavity_count,
        linker_count=linker_count,
        terminal_substituent_count=terminal_substituent_count,
        stereo_center_count=stereo_center_count,
        block_count=len(block_graph.blocks),
        has_hierarchical_blocks=has_hierarchical_blocks,
        has_block_layout_problem=has_block_layout_problem,
        preserve_only=preserve_only,
        preserve_only_required=preserve_only,
        local_repair_allowed=local_repair_allowed,
        global_redraw_allowed=global_redraw_allowed,
        internal_route_allowed=internal_route_allowed,
        reason=reason,
        block_counts=block_counts,
    )


def _normalize_atom_ids(graph: MolGraph, atom_ids: Iterable[int] | None) -> set[int]:
    if atom_ids is None:
        return set(graph.atoms)
    return {int(atom_id) for atom_id in atom_ids if int(atom_id) in graph.atoms}


def _fallback_linker_count(graph: MolGraph, selected: set[int], block_graph: object) -> int:
    rigid_blocks = [
        block
        for block in getattr(block_graph, "blocks", ()) or ()
        if getattr(block, "kind", None)
        in {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
    ]
    if len(rigid_blocks) < 2:
        return 0
    atom_to_rigid: dict[int, set[int]] = {}
    rigid_atoms: set[int] = set()
    for block in rigid_blocks:
        block_id = int(getattr(block, "id", 0) or 0)
        atoms = set(getattr(block, "atom_ids", ()) or ())
        rigid_atoms.update(atoms)
        for atom_id in atoms:
            atom_to_rigid.setdefault(atom_id, set()).add(block_id)
    linker_atoms = selected - rigid_atoms
    if not linker_atoms:
        return 0

    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in selected}
    for bond in graph.bonds.values():
        if bond.a1_id in selected and bond.a2_id in selected and bond_affects_valence(bond):
            adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
            adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)

    count = 0
    seen: set[int] = set()
    for start in sorted(linker_atoms):
        if start in seen:
            continue
        stack = [start]
        component: set[int] = set()
        adjacent_rigid: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            for neighbor in adjacency.get(atom_id, set()):
                if neighbor in linker_atoms:
                    stack.append(neighbor)
                else:
                    adjacent_rigid.update(atom_to_rigid.get(neighbor, set()))
        seen.update(component)
        if len(adjacent_rigid) >= 2:
            count += 1
    return count


def _has_intramolecular_block_layout_problem(block_graph: object, quality_report: object | None) -> bool:
    if _has_hierarchical_block_layout_signals(block_graph):
        return True
    counts = _block_kind_counts(block_graph)
    if any(
        counts.get(kind, 0) > 0
        for kind in (
            BlockKind.MACROCYCLE,
            BlockKind.CYCLOPHANE,
            BlockKind.FUSED_SYSTEM,
            BlockKind.INTRAMOLECULAR_BRIDGE,
            BlockKind.INTERNAL_CAVITY,
        )
    ):
        return True
    if counts.get(BlockKind.AROMATIC_RING, 0) >= 3:
        return True
    if counts.get(BlockKind.STEREO_CENTER, 0) > 0 and counts.get(BlockKind.AROMATIC_RING, 0) > 0:
        return True
    quality_class = getattr(quality_report, "quality_class", "") if quality_report is not None else ""
    return quality_class != "good" and bool(quality_class) and (
        counts.get(BlockKind.AROMATIC_RING, 0) >= 2
        or counts.get(BlockKind.LINKER, 0) > 0
        or counts.get(BlockKind.TERMINAL_SUBSTITUENT, 0) > 0
    )


def _has_hierarchical_block_layout_signals(block_graph: object) -> bool:
    counts = _block_kind_counts(block_graph)
    if any(
        counts.get(kind, 0) > 0
        for kind in (
            BlockKind.MACROCYCLE,
            BlockKind.CYCLOPHANE,
            BlockKind.INTERNAL_CAVITY,
            BlockKind.INTRAMOLECULAR_BRIDGE,
            BlockKind.FUSED_SYSTEM,
        )
    ):
        return True
    if counts.get(BlockKind.AROMATIC_RING, 0) >= 3:
        return True
    return counts.get(BlockKind.STEREO_CENTER, 0) >= 2


def _block_kind_counts(block_graph: object) -> dict[BlockKind, int]:
    counts: dict[BlockKind, int] = {}
    for block in getattr(block_graph, "blocks", ()) or ():
        kind = getattr(block, "kind", None)
        if isinstance(kind, BlockKind):
            counts[kind] = counts.get(kind, 0) + 1
    return counts


def _preserve_reason(**values: object) -> tuple[bool, str]:
    atom_count = int(values["atom_count"])
    ring_count = int(values["ring_count"])
    aromatic_ring_count = int(values["aromatic_ring_count"])
    stereo_center_count = int(values["stereo_center_count"])
    complex_enough_for_global_preserve = (
        atom_count >= 25
        or ring_count >= 3
        or aromatic_ring_count >= 2
        or stereo_center_count >= 2
    )
    checks = (
        (int(values["macrocycle_count"]) > 0 and complex_enough_for_global_preserve, "macrocycle"),
        (int(values["cyclophane_count"]) > 0 and complex_enough_for_global_preserve, "cyclophane"),
        (int(values["intramolecular_bridge_count"]) > 0 and complex_enough_for_global_preserve, "intramolecular_bridge"),
        (int(values["internal_cavity_count"]) > 0 and complex_enough_for_global_preserve, "internal_cavity"),
        (
            int(values["fused_system_count"]) > 0 and aromatic_ring_count >= 2,
            "fused_aromatic_systems",
        ),
        (aromatic_ring_count >= 3, "many_aromatic_rings"),
        (stereo_center_count >= 2, "multiple_stereo_centers"),
        (
            atom_count >= 45 and ring_count >= 3,
            "large_polycyclic_structure",
        ),
        (bool(values["has_hierarchical_blocks"]) and complex_enough_for_global_preserve, "hierarchical_blocks"),
        (
            bool(values["has_block_layout_problem"]) and atom_count >= 25,
            "block_layout_problem",
        ),
    )
    for matched, reason in checks:
        if matched:
            return True, reason
    return False, ""
