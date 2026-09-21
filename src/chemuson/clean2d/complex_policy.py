from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

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
            "weight": edge.weight,
            "metadata": _stable_json_value(edge.metadata),
        }
        for edge in sorted(layer_model.block_graph.edges, key=lambda item: item.id)
    ]
    block_kind_counts: dict[str, int] = {}
    for block in blocks:
        kind = str(block["kind"])
        block_kind_counts[kind] = block_kind_counts.get(kind, 0) + 1
    block_kind_counts = dict(sorted(block_kind_counts.items()))
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
    if value is None or isinstance(value, (bool, int, float, str)):
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
