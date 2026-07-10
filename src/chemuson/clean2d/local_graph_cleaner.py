from __future__ import annotations

"""Local graph-walk Clean2D cleaner.

This module intentionally avoids global redraws.  It detects local geometric
defects and accepts only small coordinate-only repairs that improve the local
score without changing the overall projection or explicit stereochemical layout.
"""

from dataclasses import dataclass, field
from enum import Enum
import math
from typing import Any, Iterable

from chemuson.clean2d.safety import ring_degeneracy_score
from chemuson.core.model import Bond, BondStyle, MolGraph, bond_affects_valence, normalize_bond_style


class LocalClean2DMode(str, Enum):
    QUICK = "quick"
    PUBLICATION = "publication"


@dataclass(frozen=True)
class LocalDefect:
    defect_type: str
    atom_ids: tuple[int, ...] = ()
    bond_ids: tuple[int, ...] = ()
    score: float = 0.0
    severity: str = "low"
    region_atom_ids: tuple[int, ...] = ()
    involves_stereo: bool = False
    involves_ring: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class LocalRepair:
    defect: LocalDefect
    coords: dict[int, tuple[float, float]]
    moved_atom_ids: tuple[int, ...]
    repair_type: str
    reason: str = ""


@dataclass(frozen=True)
class RepairValidation:
    accepted: bool
    reason: str = ""
    bond_integrity_regressions: tuple[int, ...] = ()


@dataclass(frozen=True)
class LocalClean2DReport:
    number_of_defects_detected: int = 0
    number_of_repairs_attempted: int = 0
    number_of_repairs_accepted: int = 0
    defects_by_type: dict[str, int] = field(default_factory=dict)
    rejected_repairs_by_reason: dict[str, int] = field(default_factory=dict)
    moved_atom_count: int = 0
    mean_displacement: float = 0.0
    max_displacement: float = 0.0
    bond_integrity_regressions: int = 0
    message: str = ""
    reason: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "number_of_defects_detected": self.number_of_defects_detected,
            "number_of_repairs_attempted": self.number_of_repairs_attempted,
            "number_of_repairs_accepted": self.number_of_repairs_accepted,
            "defects_by_type": dict(self.defects_by_type),
            "rejected_repairs_by_reason": dict(self.rejected_repairs_by_reason),
            "moved_atom_count": self.moved_atom_count,
            "mean_displacement": self.mean_displacement,
            "max_displacement": self.max_displacement,
            "bond_integrity_regressions": self.bond_integrity_regressions,
            "message": self.message,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class LocalClean2DResult:
    coords: dict[int, tuple[float, float]]
    changed_coords: dict[int, tuple[float, float]]
    defects: tuple[LocalDefect, ...]
    report: LocalClean2DReport

    @property
    def ok(self) -> bool:
        return self.report.number_of_repairs_accepted > 0 and bool(self.changed_coords)


@dataclass(frozen=True)
class _LocalParams:
    mode: LocalClean2DMode
    bond_low: float
    bond_high: float
    angle_error: float
    ring_score_threshold: float
    collision_distance: float
    atom_bond_distance: float
    max_passes: int
    max_repairs: int
    max_atom_displacement_total: float
    max_single_step: float
    bbox_min: float
    bbox_max: float


@dataclass(frozen=True)
class _GraphContext:
    graph: MolGraph
    atom_ids: set[int]
    bonds: list[Bond]
    coords: dict[int, tuple[float, float]]
    adjacency: dict[int, set[int]]
    bond_by_id: dict[int, Bond]
    bond_by_pair: dict[tuple[int, int], Bond]
    cycles: list[list[int]]
    aromatic_rings: list[list[int]]
    small_rings: list[list[int]]
    macrocycles: list[list[int]]
    ring_atoms: set[int]
    aromatic_atoms: set[int]
    macrocycle_atoms: set[int]
    stereo_centers: set[int]
    wedge_hash_bonds: set[int]
    rigid_scores: dict[int, float]
    rigid_blocks: tuple[frozenset[int], ...]
    atom_to_rigid_block: dict[int, frozenset[int]]


def local_graph_clean2d(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    target_bond_length: float = 42.0,
    mode: LocalClean2DMode | str = LocalClean2DMode.QUICK,
) -> LocalClean2DResult:
    """Clean a molecular graph using local, incremental graph repairs only."""
    selected = _normalize_atom_ids(graph, atom_ids)
    before = _coords_for_atoms(graph, selected)
    target = max(8.0, float(target_bond_length or 42.0))
    params = _params_for_mode(mode, target)
    if not selected or not before:
        report = LocalClean2DReport(
            message="Limpieza 2D omitida: no se detectaron defectos locales seguros",
            reason="sin_atomos",
        )
        return LocalClean2DResult(coords=before, changed_coords={}, defects=(), report=report)

    coords = dict(before)
    initial_ctx = _build_context(graph, selected, coords)
    if not initial_ctx.bonds:
        report = LocalClean2DReport(
            message="Limpieza 2D omitida: no se detectaron defectos locales seguros",
            reason="sin_enlaces",
        )
        return LocalClean2DResult(coords=coords, changed_coords={}, defects=(), report=report)

    initial_signature = stereo_layout_signature(graph, coords, selected)
    initial_crossings = _count_crossings(coords, initial_ctx.bonds)
    initial_good_aromatic_scores = _aromatic_ring_scores(initial_ctx, coords)
    initial_good_aromatic_scores = {
        key: value
        for key, value in initial_good_aromatic_scores.items()
        if value >= 0.65
    }
    initial_defects = _detect_defects(initial_ctx, target, params)
    defects_by_type = _count_defects_by_type(initial_defects)
    rejected: dict[str, int] = {}
    attempted = 0
    accepted = 0
    bond_integrity_regressions = 0

    if not initial_defects:
        report = _build_report(
            before,
            coords,
            selected,
            tuple(initial_defects),
            defects_by_type,
            rejected,
            attempted,
            accepted,
            bond_integrity_regressions=0,
            message="Limpieza 2D omitida: no se detectaron defectos locales seguros",
            reason="no_defects",
        )
        return LocalClean2DResult(coords=coords, changed_coords={}, defects=tuple(initial_defects), report=report)

    stop = False
    for _pass_idx in range(params.max_passes):
        if stop:
            break
        ctx = _build_context(graph, selected, coords)
        defects = _detect_defects(ctx, target, params)
        if not defects:
            break
        repaired_in_pass = False
        for defect in defects:
            if accepted >= params.max_repairs:
                stop = True
                break
            repairs = _repairs_for_defect(ctx, defect, coords, before, target, params)
            if not repairs:
                _inc(rejected, "sin_reparacion_local_segura")
                continue
            repair_accepted = False
            for repair in repairs:
                attempted += 1
                reason = _validate_repair(
                    graph,
                    selected,
                    ctx,
                    before,
                    coords,
                    repair.coords,
                    defect,
                    target,
                    params,
                    initial_signature,
                    initial_crossings,
                    initial_good_aromatic_scores,
                )
                if reason:
                    _inc(rejected, reason)
                    continue
                validation = validate_local_repair(
                    current=coords,
                    after=repair.coords,
                    repair=repair,
                    graph=graph,
                    atom_ids=selected,
                    bonds=ctx.bonds,
                    target=target,
                    mode=params.mode.value,
                )
                if not validation.accepted:
                    _inc(rejected, validation.reason or "reparacion_local_rechazada")
                    continue
                bond_integrity_regressions = max(
                    bond_integrity_regressions,
                    len(validation.bond_integrity_regressions),
                )
                coords = repair.coords
                accepted += 1
                repaired_in_pass = True
                repair_accepted = True
                break
            if not repair_accepted and repairs:
                continue
        if not repaired_in_pass:
            break

    changed = {
        atom_id: coords[atom_id]
        for atom_id in selected
        if atom_id in before and atom_id in coords and _distance(before[atom_id], coords[atom_id]) > 0.5
    }
    if accepted > 0 and changed:
        remaining = _detect_defects(_build_context(graph, selected, coords), target, params)
        if remaining:
            message = "Estructura 2D limpiada parcialmente"
            reason = "partial"
        else:
            message = "Estructura 2D limpiada localmente"
            reason = "cleaned"
    elif rejected.get("cambia_estereoquimica", 0) > 0:
        message = "Limpieza 2D omitida: corrección local cambiaría estereoquímica"
        reason = "stereo_constraint"
    else:
        message = "Limpieza 2D omitida: no se detectaron defectos locales seguros"
        reason = "no_safe_local_repair"

    report = _build_report(
        before,
        coords,
        selected,
        tuple(initial_defects),
        defects_by_type,
        rejected,
        attempted,
        accepted,
        bond_integrity_regressions=bond_integrity_regressions,
        message=message,
        reason=reason,
    )
    return LocalClean2DResult(coords=coords, changed_coords=changed, defects=tuple(initial_defects), report=report)


def is_complex_clean2d_graph(graph: MolGraph, atom_ids: Iterable[int] | None = None) -> bool:
    selected = _normalize_atom_ids(graph, atom_ids)
    if not selected:
        return False
    ctx = _build_context(graph, selected, _coords_for_atoms(graph, selected))
    aromatic_ring_count = len(ctx.aromatic_rings)
    if len(selected) > 25:
        return True
    if len(ctx.cycles) >= 3:
        return True
    if ctx.macrocycles:
        return True
    if aromatic_ring_count >= 2 and len(selected) > 18:
        return True
    if ctx.wedge_hash_bonds:
        return True
    return len(selected) > 18 and _has_bridged_or_fused_system(ctx.cycles)


def stereo_layout_signature(
    graph: MolGraph,
    coords: dict[int, tuple[float, float]],
    atom_ids: Iterable[int] | None = None,
) -> tuple[Any, ...]:
    selected = _normalize_atom_ids(graph, atom_ids)
    bonds = _selected_structural_bonds(graph, selected)
    adjacency = _adjacency_for_bonds(selected, bonds)
    stereo_centers = _stereo_centers(graph, selected, bonds)
    entries: list[Any] = []
    for center in sorted(stereo_centers):
        if center not in coords:
            continue
        neighbors = sorted(neigh for neigh in adjacency.get(center, set()) if neigh in coords)
        pair_signs: list[tuple[int, int, int]] = []
        cx, cy = coords[center]
        for idx, left in enumerate(neighbors):
            lx, ly = coords[left][0] - cx, coords[left][1] - cy
            llen = math.hypot(lx, ly)
            if llen <= 1e-9:
                continue
            for right in neighbors[idx + 1 :]:
                rx, ry = coords[right][0] - cx, coords[right][1] - cy
                rlen = math.hypot(rx, ry)
                if rlen <= 1e-9:
                    continue
                signed = (lx * ry - ly * rx) / (llen * rlen)
                sign = 1 if signed > 0.05 else -1 if signed < -0.05 else 0
                pair_signs.append((left, right, sign))
        entries.append(("center", center, tuple(pair_signs)))
    for bond in sorted(bonds, key=lambda item: item.id):
        style = normalize_bond_style(getattr(bond, "style", BondStyle.PLAIN))
        stereo = getattr(getattr(bond, "stereo", None), "value", getattr(bond, "stereo", None))
        if style in {BondStyle.WEDGE, BondStyle.HASHED} or (stereo and str(stereo) not in {"none", "BondStereo.NONE"}):
            entries.append(("bond", int(bond.id), int(bond.a1_id), int(bond.a2_id), getattr(style, "value", str(style)), str(stereo)))
    return tuple(entries)


def _params_for_mode(mode: LocalClean2DMode | str, target: float) -> _LocalParams:
    clean_mode = _coerce_mode(mode)
    if clean_mode == LocalClean2DMode.PUBLICATION:
        return _LocalParams(
            mode=clean_mode,
            bond_low=0.78,
            bond_high=1.25,
            angle_error=42.0,
            ring_score_threshold=0.68,
            collision_distance=0.55,
            atom_bond_distance=0.28,
            max_passes=6,
            max_repairs=80,
            max_atom_displacement_total=target * 1.5,
            max_single_step=target * 0.35,
            bbox_min=0.65,
            bbox_max=1.45,
        )
    return _LocalParams(
        mode=LocalClean2DMode.QUICK,
        bond_low=0.65,
        bond_high=1.45,
        angle_error=55.0,
        ring_score_threshold=0.55,
        collision_distance=0.45,
        atom_bond_distance=0.22,
        max_passes=3,
        max_repairs=30,
        max_atom_displacement_total=target * 0.75,
        max_single_step=target * 0.20,
        bbox_min=0.80,
        bbox_max=1.25,
    )


def _coerce_mode(mode: LocalClean2DMode | str) -> LocalClean2DMode:
    if isinstance(mode, LocalClean2DMode):
        return mode
    text = str(mode or "quick").strip().lower()
    if text in {"publication", "publish"}:
        return LocalClean2DMode.PUBLICATION
    return LocalClean2DMode.QUICK


def _build_context(
    graph: MolGraph,
    atom_ids: set[int],
    coords: dict[int, tuple[float, float]],
) -> _GraphContext:
    bonds = _selected_structural_bonds(graph, atom_ids)
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    bond_by_id = {int(bond.id): bond for bond in bonds}
    bond_by_pair = {(min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id)): bond for bond in bonds}
    cycles = _cycle_basis_ordered(atom_ids, bonds, max_size=max(3, min(len(atom_ids), 32)))
    aromatic_rings = [ring for ring in cycles if _ring_is_aromatic(ring, bond_by_pair)]
    small_rings = [ring for ring in cycles if 3 <= len(ring) <= 8]
    macrocycles = [ring for ring in cycles if len(ring) > 8]
    ring_atoms = {atom_id for ring in cycles for atom_id in ring}
    aromatic_atoms = {atom_id for ring in aromatic_rings for atom_id in ring}
    macrocycle_atoms = {atom_id for ring in macrocycles for atom_id in ring}
    wedge_hash_bonds = {
        int(bond.id)
        for bond in bonds
        if normalize_bond_style(getattr(bond, "style", BondStyle.PLAIN)) in {BondStyle.WEDGE, BondStyle.HASHED}
    }
    stereo_centers = _stereo_centers(graph, atom_ids, bonds)
    rigid_scores = _rigid_scores(atom_ids, bonds, ring_atoms, aromatic_atoms, macrocycle_atoms, stereo_centers, wedge_hash_bonds)
    rigid_blocks = _rigid_blocks(atom_ids, aromatic_rings, small_rings, stereo_centers, bonds, wedge_hash_bonds)
    atom_to_rigid_block: dict[int, frozenset[int]] = {}
    for block in rigid_blocks:
        for atom_id in block:
            atom_to_rigid_block[atom_id] = block
    return _GraphContext(
        graph=graph,
        atom_ids=set(atom_ids),
        bonds=bonds,
        coords=dict(coords),
        adjacency=adjacency,
        bond_by_id=bond_by_id,
        bond_by_pair=bond_by_pair,
        cycles=cycles,
        aromatic_rings=aromatic_rings,
        small_rings=small_rings,
        macrocycles=macrocycles,
        ring_atoms=ring_atoms,
        aromatic_atoms=aromatic_atoms,
        macrocycle_atoms=macrocycle_atoms,
        stereo_centers=stereo_centers,
        wedge_hash_bonds=wedge_hash_bonds,
        rigid_scores=rigid_scores,
        rigid_blocks=rigid_blocks,
        atom_to_rigid_block=atom_to_rigid_block,
    )


def _detect_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    defects: list[LocalDefect] = []
    defects.extend(_bond_length_defects(ctx, target, params))
    defects.extend(_angle_defects(ctx, target, params))
    defects.extend(_ring_defects(ctx, target, params))
    defects.extend(_collision_defects(ctx, target, params))
    defects.extend(_crossing_defects(ctx, target))
    defects.extend(_macrocycle_defects(ctx, target, params))
    defects.sort(key=lambda item: (_defect_priority(item), -item.score, item.atom_ids, item.bond_ids))
    return defects


def _defect_priority(defect: LocalDefect) -> int:
    return {
        "crossing": 0,
        "bond_length": 1,
        "macrocycle_segment": 1,
        "collision_atom": 2,
        "collision_atom_bond": 2,
        "aromatic_ring": 3,
        "small_ring": 3,
        "angle": 4,
    }.get(defect.defect_type, 9)


def _bond_length_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    macrocycle_edges = _ring_edge_pairs(ctx.macrocycles)
    for bond in ctx.bonds:
        if bond.a1_id not in ctx.coords or bond.a2_id not in ctx.coords:
            continue
        pair = (min(bond.a1_id, bond.a2_id), max(bond.a1_id, bond.a2_id))
        if pair in macrocycle_edges:
            continue
        length = _distance(ctx.coords[bond.a1_id], ctx.coords[bond.a2_id])
        desired = _target_length_for_bond(bond, target)
        if desired <= 1e-9:
            continue
        ratio = length / desired
        if params.bond_low <= ratio <= params.bond_high:
            continue
        score = abs(ratio - 1.0)
        atoms = (int(bond.a1_id), int(bond.a2_id))
        out.append(
            LocalDefect(
                "bond_length",
                atom_ids=atoms,
                bond_ids=(int(bond.id),),
                score=score,
                severity=_severity(score, high=0.65, medium=0.35),
                region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                involves_stereo=bool(set(atoms) & ctx.stereo_centers or int(bond.id) in ctx.wedge_hash_bonds),
                involves_ring=bool(set(atoms) & ctx.ring_atoms),
                metadata={"ratio": ratio, "length": length, "desired": desired},
            )
        )
    return out


def _angle_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    small_ring_triples = _small_ring_angle_triples(ctx.small_rings, ctx.adjacency)
    for center in sorted(ctx.atom_ids):
        if center in ctx.stereo_centers:
            continue
        if center not in ctx.coords:
            continue
        neighbors = sorted(neigh for neigh in ctx.adjacency.get(center, set()) if neigh in ctx.coords)
        if len(neighbors) < 2 or len(neighbors) >= 4:
            continue
        geometry = _atom_geometry(ctx, center)
        preferred = _preferred_angle_for_center(ctx, center, neighbors, geometry)
        for idx, left in enumerate(neighbors):
            for right in neighbors[idx + 1 :]:
                if frozenset((left, center, right)) in small_ring_triples:
                    continue
                angle = _angle_between(ctx.coords[center], ctx.coords[left], ctx.coords[right])
                diff = min(abs(angle - preferred), abs(angle - (360.0 - preferred)))
                if geometry == "sp3":
                    diff = min(diff, abs(angle - 120.0))
                if diff < params.angle_error:
                    continue
                atoms = (left, center, right)
                out.append(
                    LocalDefect(
                        "angle",
                        atom_ids=atoms,
                        score=diff / 30.0,
                        severity=_severity(diff, high=75.0, medium=params.angle_error),
                        region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                        involves_stereo=bool(set(atoms) & ctx.stereo_centers),
                        involves_ring=bool(set(atoms) & ctx.ring_atoms),
                        metadata={"angle": angle, "preferred": preferred, "diff": diff},
                    )
                )
    return out


def _ring_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    for ring in ctx.aromatic_rings:
        if not all(atom_id in ctx.coords for atom_id in ring):
            continue
        score = ring_degeneracy_score(ctx.coords, set(ring))
        length_error = _ring_length_error(ctx, ring, target * 0.97)
        if score >= params.ring_score_threshold and length_error <= 0.22:
            continue
        defect_score = max(0.0, params.ring_score_threshold - score) + length_error
        out.append(
            LocalDefect(
                "aromatic_ring",
                atom_ids=tuple(ring),
                bond_ids=tuple(_ring_bond_ids(ctx, ring)),
                score=defect_score,
                severity=_severity(defect_score, high=0.65, medium=0.30),
                region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(ring), radius=1))),
                involves_stereo=bool(set(ring) & ctx.stereo_centers),
                involves_ring=True,
                metadata={"ring_score": score, "length_error": length_error},
            )
        )
    aromatic_keys = {tuple(ring) for ring in ctx.aromatic_rings}
    for ring in ctx.small_rings:
        if tuple(ring) in aromatic_keys or not all(atom_id in ctx.coords for atom_id in ring):
            continue
        score = ring_degeneracy_score(ctx.coords, set(ring))
        if score >= 0.12:
            continue
        out.append(
            LocalDefect(
                "small_ring",
                atom_ids=tuple(ring),
                bond_ids=tuple(_ring_bond_ids(ctx, ring)),
                score=0.12 - score,
                severity="high",
                region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(ring), radius=1))),
                involves_stereo=bool(set(ring) & ctx.stereo_centers),
                involves_ring=True,
                metadata={"ring_score": score},
            )
        )
    return out


def _collision_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    bonded = {(min(b.a1_id, b.a2_id), max(b.a1_id, b.a2_id)) for b in ctx.bonds}
    ids = sorted(atom_id for atom_id in ctx.atom_ids if atom_id in ctx.coords)
    atom_threshold = target * params.collision_distance
    for idx, a_id in enumerate(ids):
        for b_id in ids[idx + 1 :]:
            if (min(a_id, b_id), max(a_id, b_id)) in bonded:
                continue
            dist = _distance(ctx.coords[a_id], ctx.coords[b_id])
            if dist >= atom_threshold:
                continue
            atoms = (a_id, b_id)
            score = (atom_threshold - dist) / target
            out.append(
                LocalDefect(
                    "collision_atom",
                    atom_ids=atoms,
                    score=score,
                    severity=_severity(score, high=0.30, medium=0.12),
                    region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                    involves_stereo=bool(set(atoms) & ctx.stereo_centers),
                    involves_ring=bool(set(atoms) & ctx.ring_atoms),
                    metadata={"distance": dist, "threshold": atom_threshold},
                )
            )

    atom_bond_threshold = target * params.atom_bond_distance
    for atom_id in ids:
        ax, ay = ctx.coords[atom_id]
        for bond in ctx.bonds:
            if atom_id in {bond.a1_id, bond.a2_id}:
                continue
            if bond.a1_id not in ctx.coords or bond.a2_id not in ctx.coords:
                continue
            dist, px, py = _point_segment_distance(ax, ay, ctx.coords[bond.a1_id], ctx.coords[bond.a2_id])
            if dist >= atom_bond_threshold:
                continue
            atoms = (atom_id, bond.a1_id, bond.a2_id)
            score = (atom_bond_threshold - dist) / target
            out.append(
                LocalDefect(
                    "collision_atom_bond",
                    atom_ids=atoms,
                    bond_ids=(int(bond.id),),
                    score=score,
                    severity=_severity(score, high=0.16, medium=0.08),
                    region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                    involves_stereo=bool(set(atoms) & ctx.stereo_centers or int(bond.id) in ctx.wedge_hash_bonds),
                    involves_ring=bool(set(atoms) & ctx.ring_atoms),
                    metadata={"distance": dist, "threshold": atom_bond_threshold, "closest": (px, py)},
                )
            )
    return out


def _crossing_defects(ctx: _GraphContext, target: float) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    for idx, b1 in enumerate(ctx.bonds):
        if b1.a1_id not in ctx.coords or b1.a2_id not in ctx.coords:
            continue
        for b2 in ctx.bonds[idx + 1 :]:
            if b2.a1_id not in ctx.coords or b2.a2_id not in ctx.coords:
                continue
            if {b1.a1_id, b1.a2_id} & {b2.a1_id, b2.a2_id}:
                continue
            if not _segments_intersect(ctx.coords[b1.a1_id], ctx.coords[b1.a2_id], ctx.coords[b2.a1_id], ctx.coords[b2.a2_id]):
                continue
            atoms = (b1.a1_id, b1.a2_id, b2.a1_id, b2.a2_id)
            out.append(
                LocalDefect(
                    "crossing",
                    atom_ids=atoms,
                    bond_ids=(int(b1.id), int(b2.id)),
                    score=4.0,
                    severity="high",
                    region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                    involves_stereo=bool(set(atoms) & ctx.stereo_centers or {int(b1.id), int(b2.id)} & ctx.wedge_hash_bonds),
                    involves_ring=bool(set(atoms) & ctx.ring_atoms),
                )
            )
    return out


def _macrocycle_defects(ctx: _GraphContext, target: float, params: _LocalParams) -> list[LocalDefect]:
    out: list[LocalDefect] = []
    seen_pairs: set[tuple[int, int]] = set()
    for ring in ctx.macrocycles:
        for idx, a_id in enumerate(ring):
            b_id = ring[(idx + 1) % len(ring)]
            pair = (min(a_id, b_id), max(a_id, b_id))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            bond = ctx.bond_by_pair.get(pair)
            if bond is None or a_id not in ctx.coords or b_id not in ctx.coords:
                continue
            desired = _target_length_for_bond(bond, target)
            length = _distance(ctx.coords[a_id], ctx.coords[b_id])
            ratio = length / max(desired, 1e-9)
            if 0.45 <= ratio <= 1.85:
                continue
            atoms = (a_id, b_id)
            out.append(
                LocalDefect(
                    "macrocycle_segment",
                    atom_ids=atoms,
                    bond_ids=(int(bond.id),),
                    score=abs(ratio - 1.0) * 0.5,
                    severity=_severity(abs(ratio - 1.0), high=0.85, medium=0.45),
                    region_atom_ids=tuple(sorted(_bfs_region(ctx.adjacency, set(atoms), radius=2))),
                    involves_stereo=bool(set(atoms) & ctx.stereo_centers or int(bond.id) in ctx.wedge_hash_bonds),
                    involves_ring=True,
                    metadata={"ratio": ratio, "length": length, "desired": desired},
                )
            )
    return out


def _repairs_for_defect(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    before: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    if defect.defect_type == "bond_length":
        repair = _repair_bond_length_local(ctx, defect, coords, target, params)
        return [repair] if repair is not None else []
    if defect.defect_type == "angle":
        return _repair_angle_local(ctx, defect, coords, target, params)
    if defect.defect_type in {"collision_atom", "collision_atom_bond"}:
        return _repair_collision_local(ctx, defect, coords, target, params)
    if defect.defect_type in {"aromatic_ring", "small_ring"}:
        repair = _repair_ring_local(ctx, defect, coords, target, params)
        return [repair] if repair is not None else []
    if defect.defect_type == "macrocycle_segment":
        repair = _repair_macrocycle_local(ctx, defect, coords, target, params)
        return [repair] if repair is not None else []
    if defect.defect_type == "crossing":
        return _repair_crossing_local(ctx, defect, coords, target, params)
    return []


def _repair_bond_length_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> LocalRepair | None:
    if not defect.bond_ids:
        return None
    bond = ctx.bond_by_id.get(defect.bond_ids[0])
    if bond is None:
        return None
    side_info = _movable_side_for_bond(ctx, bond)
    if side_info is None:
        return None
    side, anchor, root = side_info
    if anchor not in coords or root not in coords:
        return None
    ax, ay = coords[anchor]
    rx, ry = coords[root]
    dx = rx - ax
    dy = ry - ay
    dist = math.hypot(dx, dy)
    if dist <= 1e-9:
        dx, dy = _deterministic_unit_vector(anchor, root)
        dist = 1.0
    desired = _target_length_for_bond(bond, target)
    new_root = (ax + dx / dist * desired, ay + dy / dist * desired)
    shift = _clamp_vector(new_root[0] - rx, new_root[1] - ry, params.max_single_step)
    if math.hypot(*shift) <= 0.05:
        return None
    out = dict(coords)
    moved = []
    for atom_id in sorted(side):
        if atom_id not in out:
            continue
        x, y = out[atom_id]
        out[atom_id] = (x + shift[0], y + shift[1])
        moved.append(atom_id)
    return LocalRepair(defect, out, tuple(moved), "repair_bond_length_local")


def _repair_angle_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    if len(defect.atom_ids) != 3:
        return []
    left, center, right = defect.atom_ids
    if center in ctx.stereo_centers:
        return []
    preferred = float(defect.metadata.get("preferred", 120.0))
    repairs: list[LocalRepair] = []
    for moving, fixed in ((right, left), (left, right)):
        side = _component_without_edge(moving, center, ctx.adjacency, ctx.atom_ids)
        if not side or center in side or fixed in side:
            continue
        if _side_rigidity(ctx, side) > max(8.0, len(side) * 4.0):
            continue
        current_angle = _angle_deg(coords[center], coords[moving])
        fixed_angle = _angle_deg(coords[center], coords[fixed])
        candidates = (fixed_angle + preferred, fixed_angle - preferred)
        desired_angle = min(candidates, key=lambda angle: abs(_signed_angle_delta(angle, current_angle)))
        delta = _signed_angle_delta(desired_angle, current_angle)
        delta = max(-18.0, min(18.0, delta))
        if abs(delta) < 1.0:
            continue
        out = _rotate_atoms(coords, side, coords[center], delta, params.max_single_step)
        moved = tuple(sorted(atom_id for atom_id in side if _distance(coords[atom_id], out[atom_id]) > 0.05))
        if moved:
            repairs.append(LocalRepair(defect, out, moved, "repair_angle_local"))
    repairs.sort(key=lambda repair: (len(repair.moved_atom_ids), sum(ctx.rigid_scores.get(atom_id, 0.0) for atom_id in repair.moved_atom_ids)))
    return repairs[:2]


def _repair_collision_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    if defect.defect_type == "collision_atom" and len(defect.atom_ids) >= 2:
        a_id, b_id = defect.atom_ids[:2]
        threshold = float(defect.metadata.get("threshold", target * params.collision_distance))
        return _atom_pair_separation_repairs(ctx, defect, coords, a_id, b_id, threshold, params)
    if defect.defect_type == "collision_atom_bond" and len(defect.atom_ids) >= 3:
        atom_id, a_id, b_id = defect.atom_ids[:3]
        closest = defect.metadata.get("closest")
        if not isinstance(closest, tuple) or len(closest) != 2:
            ax, ay = coords[atom_id]
            _dist, px, py = _point_segment_distance(ax, ay, coords[a_id], coords[b_id])
            closest = (px, py)
        threshold = float(defect.metadata.get("threshold", target * params.atom_bond_distance))
        return _move_atom_away_repairs(ctx, defect, coords, atom_id, (float(closest[0]), float(closest[1])), threshold, params)
    return []


def _repair_ring_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> LocalRepair | None:
    ring = list(defect.atom_ids)
    if len(ring) < 3 or not all(atom_id in coords for atom_id in ring):
        return None
    for block in ctx.rigid_blocks:
        if set(ring) & set(block) and set(ring) != set(block):
            return None
    if defect.defect_type == "aromatic_ring" and defect.metadata.get("ring_score", 0.0) >= 0.75:
        return None
    center = _center({atom_id: coords[atom_id] for atom_id in ring})
    edge_length = target * (0.97 if defect.defect_type == "aromatic_ring" else 1.0)
    radius = edge_length / (2.0 * math.sin(math.pi / len(ring)))
    candidates = _regular_ring_candidates(ring, coords, center, radius)
    if not candidates:
        return None
    best = min(candidates, key=lambda item: _ring_candidate_score(ctx, ring, coords, item))
    out = dict(coords)
    deltas: dict[int, tuple[float, float]] = {}
    max_delta = 0.0
    for atom_id in ring:
        dx = best[atom_id][0] - coords[atom_id][0]
        dy = best[atom_id][1] - coords[atom_id][1]
        max_delta = max(max_delta, math.hypot(dx, dy))
        deltas[atom_id] = (dx, dy)
    if max_delta <= 0.05:
        return None
    scale = min(1.0, params.max_single_step / max_delta)
    moved: set[int] = set()
    ring_set = set(ring)
    for atom_id in ring:
        dx, dy = deltas[atom_id]
        dx *= scale
        dy *= scale
        x, y = out[atom_id]
        out[atom_id] = (x + dx, y + dy)
        moved.add(atom_id)
        for neigh in sorted(ctx.adjacency.get(atom_id, set()) - ring_set):
            side = _component_without_edge(neigh, atom_id, ctx.adjacency, ctx.atom_ids)
            if not side or side & ring_set:
                continue
            if len(side) > 6 or _side_rigidity(ctx, side) > len(side) * 2.5:
                continue
            for side_atom in side:
                sx, sy = out[side_atom]
                out[side_atom] = (sx + dx, sy + dy)
                moved.add(side_atom)
    return LocalRepair(defect, out, tuple(sorted(moved)), f"repair_{defect.defect_type}_local")


def _repair_macrocycle_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> LocalRepair | None:
    if params.mode == LocalClean2DMode.QUICK:
        return None
    if not defect.bond_ids:
        return None
    bond = ctx.bond_by_id.get(defect.bond_ids[0])
    if bond is None:
        return None
    return _repair_bond_length_local(ctx, defect, coords, target, params)


def _repair_crossing_local(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    if len(defect.bond_ids) != 2:
        return []
    bonds = [ctx.bond_by_id.get(bond_id) for bond_id in defect.bond_ids]
    if bonds[0] is None or bonds[1] is None:
        return []
    repairs: list[LocalRepair] = []
    other_midpoints = []
    for bond in bonds:
        assert bond is not None
        other_midpoints.append(((coords[bond.a1_id][0] + coords[bond.a2_id][0]) * 0.5, (coords[bond.a1_id][1] + coords[bond.a2_id][1]) * 0.5))
    for bond_idx, bond in enumerate(bonds):
        assert bond is not None
        other_mid = other_midpoints[1 - bond_idx]
        for root, anchor in ((bond.a1_id, bond.a2_id), (bond.a2_id, bond.a1_id)):
            side = _component_without_edge(root, anchor, ctx.adjacency, ctx.atom_ids)
            if not side or anchor in side:
                continue
            if _side_rigidity(ctx, side) > len(side) * 3.0 or len(side) > max(8, len(ctx.atom_ids) // 3):
                continue
            root_pt = coords[root]
            vx = root_pt[0] - other_mid[0]
            vy = root_pt[1] - other_mid[1]
            if math.hypot(vx, vy) <= 1e-9:
                vx, vy = _perpendicular(coords[bond.a1_id], coords[bond.a2_id])
            shift = _clamp_vector(vx, vy, params.max_single_step)
            out = dict(coords)
            moved = []
            for atom_id in side:
                x, y = out[atom_id]
                out[atom_id] = (x + shift[0], y + shift[1])
                moved.append(atom_id)
            repairs.append(LocalRepair(defect, out, tuple(sorted(moved)), "repair_crossing_local"))
    repairs.sort(key=lambda repair: (len(repair.moved_atom_ids), sum(ctx.rigid_scores.get(atom_id, 0.0) for atom_id in repair.moved_atom_ids)))
    return repairs[:4]


def _atom_pair_separation_repairs(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    a_id: int,
    b_id: int,
    threshold: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    repairs: list[LocalRepair] = []
    repairs.extend(_move_atom_away_repairs(ctx, defect, coords, a_id, coords[b_id], threshold, params))
    repairs.extend(_move_atom_away_repairs(ctx, defect, coords, b_id, coords[a_id], threshold, params))
    repairs.sort(key=lambda repair: (len(repair.moved_atom_ids), sum(ctx.rigid_scores.get(atom_id, 0.0) for atom_id in repair.moved_atom_ids)))
    return repairs[:4]


def _move_atom_away_repairs(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    atom_id: int,
    other_point: tuple[float, float],
    threshold: float,
    params: _LocalParams,
) -> list[LocalRepair]:
    if atom_id not in coords:
        return []
    ax, ay = coords[atom_id]
    dx = ax - other_point[0]
    dy = ay - other_point[1]
    dist = math.hypot(dx, dy)
    if dist <= 1e-9:
        dx, dy = _deterministic_unit_vector(atom_id, int(sum(defect.atom_ids) if defect.atom_ids else 1))
        dist = 1.0
    push = min(params.max_single_step, max(target_min := threshold - dist + threshold * 0.08, threshold * 0.10))
    if target_min <= 0.0:
        push = min(params.max_single_step, threshold * 0.10)
    shift = _clamp_vector(dx / dist * push, dy / dist * push, params.max_single_step)
    side_options = _movable_sides_for_atom(ctx, atom_id)
    repairs: list[LocalRepair] = []
    for side in side_options:
        out = dict(coords)
        moved = []
        for side_atom in side:
            if side_atom not in out:
                continue
            x, y = out[side_atom]
            out[side_atom] = (x + shift[0], y + shift[1])
            moved.append(side_atom)
        if moved:
            repairs.append(LocalRepair(defect, out, tuple(sorted(moved)), "repair_collision_local"))
    return repairs


def validate_local_repair(
    *,
    current: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    repair: LocalRepair,
    graph: MolGraph,
    atom_ids: Iterable[int],
    bonds: list[Bond],
    target: float,
    mode: str = "quick",
) -> RepairValidation:
    selected = set(atom_ids)
    if set(after) != set(current):
        return RepairValidation(False, "coordenadas_incompletas")
    if any(not (math.isfinite(x) and math.isfinite(y)) for x, y in after.values()):
        return RepairValidation(False, "coordenadas_no_finitas")
    if stereo_layout_signature(graph, current, selected) != stereo_layout_signature(graph, after, selected):
        return RepairValidation(False, "cambia_estereoquimica")

    ctx = _build_context(graph, selected, current)
    after_ctx = _build_context(graph, selected, after)
    moved = set(repair.moved_atom_ids) or _moved_atoms_from_coords(current, after, selected)
    if not moved:
        return RepairValidation(False, "sin_movimiento")

    if _count_crossings(after, bonds) > _count_crossings(current, bonds):
        return RepairValidation(False, "crea_cruces")

    block_reason = _validate_moved_blocks(ctx, repair, moved)
    if block_reason:
        return RepairValidation(False, block_reason)

    regressions = _bond_integrity_regressions(
        current,
        after,
        bonds,
        set(repair.defect.bond_ids),
        target,
        strict=(str(mode or "quick") == "quick"),
    )
    if regressions:
        return RepairValidation(False, "reparacion_local_produciria_enlaces_peores", tuple(regressions))

    before_global = _global_local_score(ctx, current, target, _params_for_mode(mode, target))
    after_global = _global_local_score(after_ctx, after, target, _params_for_mode(mode, target))
    if after_global > before_global + 0.05:
        return RepairValidation(False, "empeora_score_global")

    return RepairValidation(True)


def _validate_repair(
    graph: MolGraph,
    atom_ids: set[int],
    ctx: _GraphContext,
    original: dict[int, tuple[float, float]],
    current: dict[int, tuple[float, float]],
    candidate: dict[int, tuple[float, float]],
    defect: LocalDefect,
    target: float,
    params: _LocalParams,
    initial_signature: tuple[Any, ...],
    initial_crossings: int,
    initial_good_aromatic_scores: dict[tuple[int, ...], float],
) -> str:
    if set(candidate) != set(current):
        return "coordenadas_incompletas"
    if any(not (math.isfinite(x) and math.isfinite(y)) for x, y in candidate.values()):
        return "coordenadas_no_finitas"
    if stereo_layout_signature(graph, candidate, atom_ids) != initial_signature:
        return "cambia_estereoquimica"
    if _count_crossings(candidate, ctx.bonds) > max(initial_crossings, _count_crossings(current, ctx.bonds)):
        return "crea_cruces"
    for ring_key, before_score in initial_good_aromatic_scores.items():
        after_score = ring_degeneracy_score(candidate, set(ring_key))
        if after_score + 0.03 < before_score:
            return "empeora_anillo_aromatico"
    for atom_id in atom_ids:
        if atom_id not in original or atom_id not in candidate:
            continue
        if _distance(original[atom_id], candidate[atom_id]) > params.max_atom_displacement_total + 1e-6:
            return "desplazamiento_excesivo"
    bbox_ratio = _bbox_ratio(original, candidate, atom_ids)
    if bbox_ratio < params.bbox_min or bbox_ratio > params.bbox_max:
        return "cambio_proyeccion_global"
    before_local = _defect_score(ctx, defect, current, target, params)
    after_ctx = _build_context(graph, atom_ids, candidate)
    after_local = _defect_score(after_ctx, defect, candidate, target, params)
    if defect.defect_type in {"aromatic_ring", "small_ring"}:
        ring = set(defect.atom_ids)
        before_ring = ring_degeneracy_score(current, ring)
        after_ring = ring_degeneracy_score(candidate, ring)
        before_len = _ring_length_error(ctx, list(defect.atom_ids), target * 0.97, coords=current)
        after_len = _ring_length_error(after_ctx, list(defect.atom_ids), target * 0.97, coords=candidate)
        if after_ring <= before_ring + 0.02:
            return "no_mejora_local"
        if after_len > max(before_len + 0.22, 0.36):
            return "empeora_longitud_anillo"
    else:
        min_improvement = max(0.001, before_local * 0.10)
        if after_local >= before_local - min_improvement:
            return "no_mejora_local"
    before_global = _global_local_score(ctx, current, target, params)
    after_global = _global_local_score(after_ctx, candidate, target, params)
    if after_global > before_global + max(0.05, before_local * 0.25):
        return "empeora_geometria_global"
    return ""


def _validate_moved_blocks(ctx: _GraphContext, repair: LocalRepair, moved: set[int]) -> str:
    if len(moved) == 1:
        atom_id = next(iter(moved))
        if atom_id in ctx.ring_atoms:
            return "movimiento_atomico_anillo"
        if atom_id in ctx.stereo_centers:
            return "movimiento_atomico_estereo"
        if _structural_degree(ctx, atom_id) > 1:
            return "movimiento_atomico_no_terminal"

    target_atoms = set(repair.defect.atom_ids)
    ring_repair = repair.repair_type in {"repair_aromatic_ring_local", "repair_small_ring_local"}
    for block in ctx.rigid_blocks:
        moved_part = moved & set(block)
        if not moved_part:
            continue
        if moved_part != set(block):
            return "movimiento_incoherente_bloque_rigido"
        if not ring_repair and not set(block) <= target_atoms:
            return "mueve_bloque_rigido_sin_razon"
    if moved & ctx.ring_atoms and not ring_repair:
        return "movimiento_anillo_no_permitido"
    return ""


def _moved_atoms_from_coords(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> set[int]:
    return {
        atom_id
        for atom_id in atom_ids
        if atom_id in before and atom_id in after and _distance(before[atom_id], after[atom_id]) > 0.5
    }


def _structural_degree(ctx: _GraphContext, atom_id: int) -> int:
    return len(ctx.adjacency.get(atom_id, set()) & ctx.atom_ids)


def _bond_integrity_regressions(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target_bond_ids: set[int],
    target: float,
    *,
    strict: bool,
) -> list[int]:
    regressions: list[int] = []
    low = 0.70
    high = 1.35 if strict else 1.50
    for bond in bonds:
        if bond.a1_id not in before or bond.a2_id not in before:
            continue
        if bond.a1_id not in after or bond.a2_id not in after:
            regressions.append(int(bond.id))
            continue
        desired = _target_length_for_bond(bond, target)
        before_len = _distance(before[bond.a1_id], before[bond.a2_id])
        after_len = _distance(after[bond.a1_id], after[bond.a2_id])
        before_ratio = before_len / max(desired, 1e-9)
        after_ratio = after_len / max(desired, 1e-9)
        before_good = low <= before_ratio <= high
        if int(bond.id) in target_bond_ids:
            if after_ratio < 0.55 or after_len >= before_len - 0.1:
                regressions.append(int(bond.id))
            continue
        if before_good and not (low <= after_ratio <= high):
            regressions.append(int(bond.id))
            continue
        if after_len > before_len * 1.20 and after_ratio > before_ratio + 0.08:
            regressions.append(int(bond.id))
    return regressions


def _defect_score(
    ctx: _GraphContext,
    defect: LocalDefect,
    coords: dict[int, tuple[float, float]],
    target: float,
    params: _LocalParams,
) -> float:
    if defect.defect_type in {"bond_length", "macrocycle_segment"} and defect.bond_ids:
        bond = ctx.bond_by_id.get(defect.bond_ids[0])
        if bond is None or bond.a1_id not in coords or bond.a2_id not in coords:
            return math.inf
        desired = _target_length_for_bond(bond, target)
        return abs(_distance(coords[bond.a1_id], coords[bond.a2_id]) - desired) / max(target, 1e-9)
    if defect.defect_type == "angle" and len(defect.atom_ids) == 3:
        left, center, right = defect.atom_ids
        if not {left, center, right} <= set(coords):
            return math.inf
        preferred = float(defect.metadata.get("preferred", 120.0))
        angle = _angle_between(coords[center], coords[left], coords[right])
        diff = min(abs(angle - preferred), abs(angle - (360.0 - preferred)))
        if _atom_geometry(ctx, center) == "sp3":
            diff = min(diff, abs(angle - 120.0))
        return diff / 30.0
    if defect.defect_type == "collision_atom" and len(defect.atom_ids) >= 2:
        a_id, b_id = defect.atom_ids[:2]
        threshold = float(defect.metadata.get("threshold", target * params.collision_distance))
        if a_id not in coords or b_id not in coords:
            return math.inf
        return max(0.0, threshold - _distance(coords[a_id], coords[b_id])) / max(target, 1e-9)
    if defect.defect_type == "collision_atom_bond" and len(defect.atom_ids) >= 3:
        atom_id, a_id, b_id = defect.atom_ids[:3]
        threshold = float(defect.metadata.get("threshold", target * params.atom_bond_distance))
        if atom_id not in coords or a_id not in coords or b_id not in coords:
            return math.inf
        dist, _px, _py = _point_segment_distance(coords[atom_id][0], coords[atom_id][1], coords[a_id], coords[b_id])
        return max(0.0, threshold - dist) / max(target, 1e-9)
    if defect.defect_type in {"aromatic_ring", "small_ring"}:
        ring = set(defect.atom_ids)
        if not ring <= set(coords):
            return math.inf
        if defect.defect_type == "aromatic_ring":
            length_error = _ring_length_error(ctx, list(defect.atom_ids), target * 0.97, coords=coords)
            return max(0.0, params.ring_score_threshold - ring_degeneracy_score(coords, ring)) + length_error
        return max(0.0, 0.12 - ring_degeneracy_score(coords, ring))
    if defect.defect_type == "crossing" and len(defect.bond_ids) == 2:
        b1 = ctx.bond_by_id.get(defect.bond_ids[0])
        b2 = ctx.bond_by_id.get(defect.bond_ids[1])
        if b1 is None or b2 is None:
            return math.inf
        if not {b1.a1_id, b1.a2_id, b2.a1_id, b2.a2_id} <= set(coords):
            return math.inf
        return 1.0 if _segments_intersect(coords[b1.a1_id], coords[b1.a2_id], coords[b2.a1_id], coords[b2.a2_id]) else 0.0
    return defect.score


def _global_local_score(ctx: _GraphContext, coords: dict[int, tuple[float, float]], target: float, params: _LocalParams) -> float:
    score = 0.0
    for bond in ctx.bonds:
        if bond.a1_id not in coords or bond.a2_id not in coords:
            score += 100.0
            continue
        desired = _target_length_for_bond(bond, target)
        ratio = _distance(coords[bond.a1_id], coords[bond.a2_id]) / max(desired, 1e-9)
        if ratio < params.bond_low or ratio > params.bond_high:
            score += abs(ratio - 1.0) ** 2 * 10.0
    atom_threshold = target * params.collision_distance
    bonded = {(min(b.a1_id, b.a2_id), max(b.a1_id, b.a2_id)) for b in ctx.bonds}
    ids = sorted(atom_id for atom_id in ctx.atom_ids if atom_id in coords)
    for idx, a_id in enumerate(ids):
        for b_id in ids[idx + 1 :]:
            if (min(a_id, b_id), max(a_id, b_id)) in bonded:
                continue
            dist = _distance(coords[a_id], coords[b_id])
            if dist < atom_threshold:
                score += ((atom_threshold - dist) / target) ** 2 * 8.0
    score += _count_crossings(coords, ctx.bonds) * 20.0
    for ring in ctx.aromatic_rings:
        score += max(0.0, params.ring_score_threshold - ring_degeneracy_score(coords, set(ring))) * 8.0
    return score


def _movable_side_for_bond(ctx: _GraphContext, bond: Bond) -> tuple[set[int], int, int] | None:
    if normalize_bond_style(getattr(bond, "style", BondStyle.PLAIN)) in {BondStyle.WEDGE, BondStyle.HASHED}:
        return None
    if int(bond.id) in ctx.wedge_hash_bonds:
        return None
    a_id, b_id = int(bond.a1_id), int(bond.a2_id)
    pair = (min(a_id, b_id), max(a_id, b_id))
    if pair in _ring_edge_pairs(ctx.small_rings):
        return None
    left = _component_without_edge(a_id, b_id, ctx.adjacency, ctx.atom_ids)
    right = _component_without_edge(b_id, a_id, ctx.adjacency, ctx.atom_ids)
    if not left or not right or left & right:
        return None
    options = [(left, b_id, a_id), (right, a_id, b_id)]
    filtered = []
    for side, anchor, root in options:
        if anchor in side or root not in side:
            continue
        if side & (ctx.ring_atoms | ctx.stereo_centers):
            continue
        if any(set(block) & side for block in ctx.rigid_blocks):
            continue
        anchor_bonus = 0.0
        side_rigidity = _side_rigidity(ctx, side)
        if side_rigidity > max(4.0, len(side) * 2.0):
            continue
        filtered.append((side_rigidity + len(side) * 0.25 + anchor_bonus, len(side), min(side), side, anchor, root))
    if not filtered:
        return None
    _score, _size, _min_id, side, anchor, root = min(filtered, key=lambda item: item[:3])
    return set(side), anchor, root


def _movable_sides_for_atom(ctx: _GraphContext, atom_id: int) -> list[set[int]]:
    sides: list[tuple[float, int, set[int]]] = []
    for neigh in sorted(ctx.adjacency.get(atom_id, set())):
        side = _component_without_edge(atom_id, neigh, ctx.adjacency, ctx.atom_ids)
        if not side or neigh in side:
            continue
        if side & (ctx.ring_atoms | ctx.stereo_centers):
            continue
        if any(set(block) & side for block in ctx.rigid_blocks):
            continue
        rigidity = _side_rigidity(ctx, side)
        if atom_id in ctx.stereo_centers and len(side) > 1:
            continue
        if rigidity > max(10.0, len(side) * 4.0):
            continue
        sides.append((rigidity + len(side) * 0.25, len(side), side))
    if (
        not sides
        and atom_id not in ctx.ring_atoms
        and atom_id not in ctx.stereo_centers
        and _structural_degree(ctx, atom_id) <= 1
    ):
        sides.append((ctx.rigid_scores.get(atom_id, 0.0), 1, {atom_id}))
    sides.sort(key=lambda item: (item[0], item[1], min(item[2])))
    return [set(side) for _score, _size, side in sides[:3]]


def _regular_ring_candidates(
    ring: list[int],
    coords: dict[int, tuple[float, float]],
    center: tuple[float, float],
    radius: float,
) -> list[dict[int, tuple[float, float]]]:
    candidates: list[dict[int, tuple[float, float]]] = []
    n = len(ring)
    for sign in (1.0, -1.0):
        sum_cos = 0.0
        sum_sin = 0.0
        for idx, atom_id in enumerate(ring):
            angle = math.atan2(coords[atom_id][1] - center[1], coords[atom_id][0] - center[0])
            base = angle - sign * 2.0 * math.pi * idx / n
            sum_cos += math.cos(base)
            sum_sin += math.sin(base)
        base_angle = math.atan2(sum_sin, sum_cos)
        candidate = {}
        for idx, atom_id in enumerate(ring):
            angle = base_angle + sign * 2.0 * math.pi * idx / n
            candidate[atom_id] = (center[0] + radius * math.cos(angle), center[1] + radius * math.sin(angle))
        candidates.append(candidate)
    return candidates


def _ring_candidate_score(
    ctx: _GraphContext,
    ring: list[int],
    coords: dict[int, tuple[float, float]],
    candidate: dict[int, tuple[float, float]],
) -> tuple[float, float]:
    rms = math.sqrt(sum(_distance(coords[atom_id], candidate[atom_id]) ** 2 for atom_id in ring) / len(ring))
    center = _center({atom_id: coords[atom_id] for atom_id in ring})
    ring_set = set(ring)
    outward_penalty = 0.0
    for atom_id in ring:
        for neigh in ctx.adjacency.get(atom_id, set()) - ring_set:
            if neigh not in coords:
                continue
            old_branch = (coords[neigh][0] - coords[atom_id][0], coords[neigh][1] - coords[atom_id][1])
            new_outward = (candidate[atom_id][0] - center[0], candidate[atom_id][1] - center[1])
            dot = old_branch[0] * new_outward[0] + old_branch[1] * new_outward[1]
            if dot < 0.0:
                outward_penalty += abs(dot) * 0.01
    return outward_penalty, rms


def _build_report(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
    defects: tuple[LocalDefect, ...],
    defects_by_type: dict[str, int],
    rejected: dict[str, int],
    attempted: int,
    accepted: int,
    bond_integrity_regressions: int,
    *,
    message: str,
    reason: str,
) -> LocalClean2DReport:
    displacements = [
        _distance(before[atom_id], after[atom_id])
        for atom_id in atom_ids
        if atom_id in before and atom_id in after
    ]
    moved = [value for value in displacements if value > 0.5]
    return LocalClean2DReport(
        number_of_defects_detected=len(defects),
        number_of_repairs_attempted=attempted,
        number_of_repairs_accepted=accepted,
        defects_by_type=dict(defects_by_type),
        rejected_repairs_by_reason=dict(rejected),
        moved_atom_count=len(moved),
        mean_displacement=(sum(displacements) / len(displacements) if displacements else 0.0),
        max_displacement=(max(displacements) if displacements else 0.0),
        bond_integrity_regressions=bond_integrity_regressions,
        message=message,
        reason=reason,
    )


def _count_defects_by_type(defects: list[LocalDefect]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for defect in defects:
        _inc(counts, defect.defect_type)
    return counts


def _inc(mapping: dict[str, int], key: str) -> None:
    mapping[key] = mapping.get(key, 0) + 1


def _normalize_atom_ids(graph: MolGraph, atom_ids: Iterable[int] | None) -> set[int]:
    if atom_ids is None:
        return set(graph.atoms.keys())
    return {int(atom_id) for atom_id in atom_ids if int(atom_id) in graph.atoms}


def _coords_for_atoms(graph: MolGraph, atom_ids: set[int]) -> dict[int, tuple[float, float]]:
    return {
        atom_id: (float(graph.atoms[atom_id].x), float(graph.atoms[atom_id].y))
        for atom_id in atom_ids
        if atom_id in graph.atoms
    }


def _selected_structural_bonds(graph: MolGraph, atom_ids: set[int]) -> list[Bond]:
    return [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in atom_ids and bond.a2_id in atom_ids and bond_affects_valence(bond)
    ]


def _adjacency_for_bonds(atom_ids: set[int], bonds: list[Bond]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for bond in bonds:
        if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
            adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
            adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    return adjacency


def _cycle_basis_ordered(atom_ids: set[int], bonds: list[Bond], *, max_size: int) -> list[list[int]]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    rings: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for bond in sorted(bonds, key=lambda item: item.id):
        path = _shortest_path_without_edge(adjacency, bond.a1_id, bond.a2_id, max_size=max_size)
        if path is None or not (3 <= len(path) <= max_size):
            continue
        key = tuple(sorted(path))
        if key in seen:
            continue
        seen.add(key)
        rings.append(path)
    rings.sort(key=lambda ring: (len(ring), min(ring), tuple(ring)))
    return rings


def _shortest_path_without_edge(
    adjacency: dict[int, set[int]],
    start: int,
    end: int,
    *,
    max_size: int,
) -> list[int] | None:
    blocked = (min(start, end), max(start, end))
    queue: list[tuple[int, list[int]]] = [(start, [start])]
    visited: set[tuple[int, int]] = set()
    while queue:
        node, path = queue.pop(0)
        if len(path) > max_size:
            continue
        if node == end and len(path) >= 3:
            return path
        state = (node, len(path))
        if state in visited:
            continue
        visited.add(state)
        for neigh in sorted(adjacency.get(node, set())):
            edge = (min(node, neigh), max(node, neigh))
            if edge == blocked or neigh in path:
                continue
            queue.append((neigh, path + [neigh]))
    return None


def _ring_is_aromatic(ring: list[int], bond_by_pair: dict[tuple[int, int], Bond]) -> bool:
    if len(ring) < 5 or len(ring) > 7:
        return False
    aromatic_edges = 0
    for idx, a_id in enumerate(ring):
        b_id = ring[(idx + 1) % len(ring)]
        bond = bond_by_pair.get((min(a_id, b_id), max(a_id, b_id)))
        if bond is not None and bool(getattr(bond, "is_aromatic", False)):
            aromatic_edges += 1
    return aromatic_edges >= max(3, len(ring) - 1)


def _ring_edge_pairs(rings: list[list[int]]) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for ring in rings:
        for idx, atom_id in enumerate(ring):
            nxt = ring[(idx + 1) % len(ring)]
            pairs.add((min(atom_id, nxt), max(atom_id, nxt)))
    return pairs


def _ring_bond_ids(ctx: _GraphContext, ring: list[int]) -> list[int]:
    ids: list[int] = []
    for idx, atom_id in enumerate(ring):
        nxt = ring[(idx + 1) % len(ring)]
        bond = ctx.bond_by_pair.get((min(atom_id, nxt), max(atom_id, nxt)))
        if bond is not None:
            ids.append(int(bond.id))
    return ids


def _stereo_centers(graph: MolGraph, atom_ids: set[int], bonds: list[Bond]) -> set[int]:
    centers: set[int] = set()
    for atom_id in atom_ids:
        atom = graph.atoms.get(atom_id)
        if atom is None:
            continue
        if any(getattr(atom, attr, None) for attr in ("stereo_cip", "stereo_axial", "stereo_helical", "stereo_si_re")):
            centers.add(atom_id)
    for bond in bonds:
        style = normalize_bond_style(getattr(bond, "style", BondStyle.PLAIN))
        stereo = getattr(getattr(bond, "stereo", None), "value", getattr(bond, "stereo", None))
        has_stereo = style in {BondStyle.WEDGE, BondStyle.HASHED} or (
            stereo and str(stereo) not in {"none", "BondStereo.NONE"}
        )
        has_stereo = has_stereo or any(
            getattr(bond, attr, None)
            for attr in ("stereo_axial", "stereo_endo_exo", "stereo_helical")
        )
        if has_stereo:
            centers.update({int(bond.a1_id), int(bond.a2_id)})
    return centers & atom_ids


def _rigid_scores(
    atom_ids: set[int],
    bonds: list[Bond],
    ring_atoms: set[int],
    aromatic_atoms: set[int],
    macrocycle_atoms: set[int],
    stereo_centers: set[int],
    wedge_hash_bonds: set[int],
) -> dict[int, float]:
    scores = {atom_id: 0.0 for atom_id in atom_ids}
    for atom_id in ring_atoms:
        scores[atom_id] = scores.get(atom_id, 0.0) + 4.0
    for atom_id in aromatic_atoms:
        scores[atom_id] = scores.get(atom_id, 0.0) + 6.0
    for atom_id in macrocycle_atoms:
        scores[atom_id] = scores.get(atom_id, 0.0) + 1.5
    for atom_id in stereo_centers:
        scores[atom_id] = scores.get(atom_id, 0.0) + 8.0
    for bond in bonds:
        if int(bond.id) in wedge_hash_bonds:
            scores[bond.a1_id] = scores.get(bond.a1_id, 0.0) + 5.0
            scores[bond.a2_id] = scores.get(bond.a2_id, 0.0) + 5.0
    return scores


def _rigid_blocks(
    atom_ids: set[int],
    aromatic_rings: list[list[int]],
    small_rings: list[list[int]],
    stereo_centers: set[int],
    bonds: list[Bond],
    wedge_hash_bonds: set[int],
) -> tuple[frozenset[int], ...]:
    blocks: list[set[int]] = []
    for ring in aromatic_rings:
        _merge_block(blocks, set(ring), min_shared=2)
    aromatic_atoms = {atom_id for ring in aromatic_rings for atom_id in ring}
    for ring in small_rings:
        ring_set = set(ring)
        if ring_set <= aromatic_atoms:
            continue
        _merge_block(blocks, ring_set, min_shared=1)
    for center in stereo_centers:
        block = {center}
        for bond in bonds:
            if center in {bond.a1_id, bond.a2_id} and int(bond.id) in wedge_hash_bonds:
                block.update({bond.a1_id, bond.a2_id})
        _merge_block(blocks, block & atom_ids, min_shared=1)
    return tuple(frozenset(block) for block in blocks if block)


def _merge_block(blocks: list[set[int]], new_block: set[int], *, min_shared: int) -> None:
    if not new_block:
        return
    merged = set(new_block)
    keep: list[set[int]] = []
    for block in blocks:
        if len(block & merged) >= min_shared:
            merged.update(block)
        else:
            keep.append(block)
    keep.append(merged)
    blocks[:] = keep


def _has_bridged_or_fused_system(cycles: list[list[int]]) -> bool:
    if len(cycles) < 2:
        return False
    for idx, left in enumerate(cycles):
        left_set = set(left)
        for right in cycles[idx + 1 :]:
            shared = left_set & set(right)
            if len(shared) >= 2:
                return True
    return False


def _bfs_region(adjacency: dict[int, set[int]], starts: set[int], radius: int) -> set[int]:
    region = set(starts)
    frontier = set(starts)
    for _ in range(radius):
        next_frontier: set[int] = set()
        for atom_id in frontier:
            next_frontier.update(adjacency.get(atom_id, set()))
        next_frontier -= region
        if not next_frontier:
            break
        region.update(next_frontier)
        frontier = next_frontier
    return region


def _component_without_edge(
    start: int,
    blocked_neighbor: int,
    adjacency: dict[int, set[int]],
    atom_ids: set[int],
) -> set[int]:
    visited: set[int] = set()
    stack = [start]
    blocked = (min(start, blocked_neighbor), max(start, blocked_neighbor))
    while stack:
        node = stack.pop()
        if node in visited or node not in atom_ids:
            continue
        visited.add(node)
        for neigh in adjacency.get(node, set()):
            edge = (min(node, neigh), max(node, neigh))
            if edge == blocked or neigh in visited:
                continue
            stack.append(neigh)
    return visited


def _side_rigidity(ctx: _GraphContext, side: set[int]) -> float:
    return sum(ctx.rigid_scores.get(atom_id, 0.0) for atom_id in side)


def _small_ring_angle_triples(rings: list[list[int]], adjacency: dict[int, set[int]]) -> set[frozenset[int]]:
    triples: set[frozenset[int]] = set()
    for ring in rings:
        if len(ring) > 4:
            continue
        ring_set = set(ring)
        for center in ring:
            neighbors = sorted(adjacency.get(center, set()) & ring_set)
            for idx, left in enumerate(neighbors):
                for right in neighbors[idx + 1 :]:
                    triples.add(frozenset((left, center, right)))
    return triples


def _atom_geometry(ctx: _GraphContext, atom_id: int) -> str:
    has_triple = False
    has_double = False
    degree = 0
    for bond in ctx.bonds:
        if atom_id not in {bond.a1_id, bond.a2_id}:
            continue
        degree += 1
        order = int(getattr(bond, "order", 1) or 1)
        if order >= 3:
            has_triple = True
        if order == 2 or bool(getattr(bond, "is_aromatic", False)):
            has_double = True
    if has_triple:
        return "sp"
    if has_double:
        return "sp2"
    atom = ctx.graph.atoms.get(atom_id)
    if atom is not None and atom.element in {"B", "N", "O"} and degree <= 2:
        return "sp2"
    return "sp3"


def _preferred_angle_for_center(ctx: _GraphContext, center: int, neighbors: list[int], geometry: str) -> float:
    preferred = 180.0 if geometry == "sp" else 120.0
    if geometry == "sp3":
        preferred = 109.5
    if len(neighbors) == 2:
        for neigh in neighbors:
            bond = ctx.bond_by_pair.get((min(center, neigh), max(center, neigh)))
            if bond is not None and bool(getattr(bond, "is_aromatic", False)):
                return 120.0
    return preferred


def _target_length_for_bond(bond: Bond | None, target: float) -> float:
    if bond is None:
        return target
    order = int(getattr(bond, "order", 1) or 1)
    if order >= 3:
        return target * 0.92
    if order == 2 or bool(getattr(bond, "is_aromatic", False)):
        return target * 0.97
    return target


def _ring_length_error(
    ctx: _GraphContext,
    ring: list[int],
    desired: float,
    *,
    coords: dict[int, tuple[float, float]] | None = None,
) -> float:
    positions = coords or ctx.coords
    errors = []
    for idx, atom_id in enumerate(ring):
        nxt = ring[(idx + 1) % len(ring)]
        if atom_id not in positions or nxt not in positions:
            continue
        errors.append(abs(_distance(positions[atom_id], positions[nxt]) - desired) / max(desired, 1e-9))
    return max(errors) if errors else 0.0


def _aromatic_ring_scores(ctx: _GraphContext, coords: dict[int, tuple[float, float]]) -> dict[tuple[int, ...], float]:
    return {
        tuple(ring): ring_degeneracy_score(coords, set(ring))
        for ring in ctx.aromatic_rings
        if all(atom_id in coords for atom_id in ring)
    }


def _rotate_atoms(
    coords: dict[int, tuple[float, float]],
    atom_ids: set[int],
    origin: tuple[float, float],
    delta_deg: float,
    max_step: float,
) -> dict[int, tuple[float, float]]:
    rad = math.radians(delta_deg)
    cos_t = math.cos(rad)
    sin_t = math.sin(rad)
    raw: dict[int, tuple[float, float]] = {}
    max_delta = 0.0
    for atom_id in atom_ids:
        x, y = coords[atom_id]
        dx = x - origin[0]
        dy = y - origin[1]
        nx = origin[0] + dx * cos_t - dy * sin_t
        ny = origin[1] + dx * sin_t + dy * cos_t
        raw[atom_id] = (nx, ny)
        max_delta = max(max_delta, math.hypot(nx - x, ny - y))
    scale = min(1.0, max_step / max_delta) if max_delta > 1e-9 else 1.0
    out = dict(coords)
    for atom_id, (nx, ny) in raw.items():
        x, y = coords[atom_id]
        out[atom_id] = (x + (nx - x) * scale, y + (ny - y) * scale)
    return out


def _clamp_vector(dx: float, dy: float, max_len: float) -> tuple[float, float]:
    length = math.hypot(dx, dy)
    if length <= max_len or length <= 1e-9:
        return dx, dy
    scale = max_len / length
    return dx * scale, dy * scale


def _deterministic_unit_vector(a_id: int, b_id: int) -> tuple[float, float]:
    angle = math.radians(((int(a_id) * 37 + int(b_id) * 17) % 360) or 13)
    return math.cos(angle), math.sin(angle)


def _perpendicular(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    length = math.hypot(dx, dy)
    if length <= 1e-9:
        return 1.0, 0.0
    return -dy / length, dx / length


def _bbox_ratio(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> float:
    before_diag = _bbox_diag(before, atom_ids)
    after_diag = _bbox_diag(after, atom_ids)
    if before_diag <= 1e-9:
        return 1.0
    return after_diag / before_diag


def _bbox_diag(coords: dict[int, tuple[float, float]], atom_ids: set[int]) -> float:
    xs = [coords[atom_id][0] for atom_id in atom_ids if atom_id in coords]
    ys = [coords[atom_id][1] for atom_id in atom_ids if atom_id in coords]
    if not xs or not ys:
        return 1.0
    return math.hypot((max(xs) - min(xs)) or 1.0, (max(ys) - min(ys)) or 1.0)


def _center(coords: dict[Any, tuple[float, float]]) -> tuple[float, float]:
    if not coords:
        return 0.0, 0.0
    return (
        sum(float(x) for x, _y in coords.values()) / len(coords),
        sum(float(y) for _x, y in coords.values()) / len(coords),
    )


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(float(a[0]) - float(b[0]), float(a[1]) - float(b[1]))


def _angle_deg(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 360.0


def _angle_between(center: tuple[float, float], left: tuple[float, float], right: tuple[float, float]) -> float:
    return abs((_angle_deg(center, right) - _angle_deg(center, left) + 180.0) % 360.0 - 180.0)


def _signed_angle_delta(target: float, current: float) -> float:
    return (target - current + 180.0) % 360.0 - 180.0


def _point_segment_distance(
    px: float,
    py: float,
    a: tuple[float, float],
    b: tuple[float, float],
) -> tuple[float, float, float]:
    ax, ay = a
    bx, by = b
    dx = bx - ax
    dy = by - ay
    denom = dx * dx + dy * dy
    if denom <= 1e-12:
        return math.hypot(px - ax, py - ay), ax, ay
    t = ((px - ax) * dx + (py - ay) * dy) / denom
    t = max(0.0, min(1.0, t))
    cx = ax + t * dx
    cy = ay + t * dy
    return math.hypot(px - cx, py - cy), cx, cy


def _count_crossings(coords: dict[int, tuple[float, float]], bonds: list[Bond]) -> int:
    crossings = 0
    for idx, b1 in enumerate(bonds):
        if b1.a1_id not in coords or b1.a2_id not in coords:
            continue
        for b2 in bonds[idx + 1 :]:
            if b2.a1_id not in coords or b2.a2_id not in coords:
                continue
            if {b1.a1_id, b1.a2_id} & {b2.a1_id, b2.a2_id}:
                continue
            if _segments_intersect(coords[b1.a1_id], coords[b1.a2_id], coords[b2.a1_id], coords[b2.a2_id]):
                crossings += 1
    return crossings


def _segments_intersect(
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    p4: tuple[float, float],
) -> bool:
    def orient(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> bool:
        return min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])

    o1 = orient(p1, p2, p3)
    o2 = orient(p1, p2, p4)
    o3 = orient(p3, p4, p1)
    o4 = orient(p3, p4, p2)
    eps = 1e-9
    if abs(o1) <= eps and on_segment(p1, p2, p3):
        return False
    if abs(o2) <= eps and on_segment(p1, p2, p4):
        return False
    if abs(o3) <= eps and on_segment(p3, p4, p1):
        return False
    if abs(o4) <= eps and on_segment(p3, p4, p2):
        return False
    return (o1 > 0.0) != (o2 > 0.0) and (o3 > 0.0) != (o4 > 0.0)


def _severity(value: float, *, high: float, medium: float) -> str:
    if value >= high:
        return "high"
    if value >= medium:
        return "medium"
    return "low"
