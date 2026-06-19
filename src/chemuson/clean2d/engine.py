from __future__ import annotations

"""Unified Clean2D depiction engine.

The engine is intentionally pure: it never mutates a ``MolGraph``.  Backends
produce coordinate candidates, the ranking stage validates invariants and
visual quality, and the GUI applies the chosen coordinates as one undoable
command.
"""

from dataclasses import dataclass, field
from enum import Enum
import hashlib
import math
from typing import Any, Iterable

from chemuson.clean2d.length_only import (
    length_only_polish,
    structure_preserving_geometry_polish,
    structure_preserving_length_polish,
)
from chemuson.clean2d.local_graph_cleaner import (
    LocalClean2DMode,
    is_complex_clean2d_graph,
    local_graph_clean2d,
    stereo_layout_signature,
)
from chemuson.clean2d.safety import (
    Clean2DQualityReport,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    min_nonbonded_distance,
    ring_degeneracy_score,
)
from chemuson.clean2d.v2 import Clean2DParameters, optimize_clean2d_positions
from chemuson.core.layers import BlockEdgeKind, BlockKind, LayoutConstraint, LayoutConstraintKind, build_multilayer_chemical_graph
from chemuson.core.model import Bond, MolGraph, bond_affects_valence


class Clean2DMode(str, Enum):
    QUICK = "quick"
    PUBLICATION = "publication"
    PROPOSE = "propose"


class Clean2DInvariantError(AssertionError):
    """Raised when a candidate changes anything except coordinates."""


@dataclass(frozen=True)
class Clean2DCandidate:
    source: str
    coords: dict[int, tuple[float, float]]
    message: str = ""
    score: float = math.inf
    novelty: float = 0.0
    report: Clean2DQualityReport | None = None
    rejected: bool = False
    rejection_reason: str = ""
    geometry_hash: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Clean2DResult:
    mode: Clean2DMode
    atom_ids: set[int]
    before_coords: dict[int, tuple[float, float]]
    candidates: tuple[Clean2DCandidate, ...]
    selected: Clean2DCandidate | None
    rejected: tuple[Clean2DCandidate, ...] = ()
    message: str = ""
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.selected is not None and not self.selected.rejected


@dataclass(frozen=True)
class Clean2DLayoutQualityReport:
    quality_class: str
    reason: str = ""
    length_rms_error: float = 0.0
    length_max_error: float = 0.0
    angle_rms_deviation: float = 0.0
    angle_max_deviation: float = 0.0
    angle_penalty: float = 0.0
    severe_angle_count: int = 0
    crossings: int = 0
    min_nonbonded_distance: float = math.inf
    min_ring_degeneracy: float = math.inf
    exocyclic_ring_angle_min: float = math.inf
    bad_exocyclic_ring_count: int = 0
    visual_score: float = 0.0


@dataclass(frozen=True)
class Clean2DGraphSnapshot:
    atom_ids: frozenset[int]
    bond_ids: frozenset[int]
    atom_data: tuple[tuple[int, tuple[Any, ...]], ...]
    bond_data: tuple[tuple[int, tuple[Any, ...]], ...]
    components: tuple[tuple[int, ...], ...]


def run_clean2d_engine(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    mode: Clean2DMode | str = Clean2DMode.QUICK,
    target_bond_length: float = 42.0,
    avoid_hashes: set[str] | None = None,
    seed: int | None = None,
    rdkit_timeout_s: float = 8.0,
) -> Clean2DResult:
    mode = _coerce_mode(mode)
    selected = _normalize_atom_ids(graph, atom_ids)
    before = _coords_for_atoms(graph, selected)
    target = max(8.0, float(target_bond_length or 42.0))
    layer_model = build_multilayer_chemical_graph(graph, selected)
    has_interaction_constraints = _interaction_constraint_count(layer_model.layout_constraint_graph.constraints) > 0
    quality = classify_clean2d_layout_quality(
        graph,
        selected,
        coords=before,
        target_bond_length=target,
    )
    has_block_layout_problem = has_intramolecular_block_layout_problem(layer_model.block_graph, quality)
    has_hierarchical_blocks = has_hierarchical_block_layout_signals(layer_model.block_graph)
    if (
        mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}
        and not has_interaction_constraints
        and not has_block_layout_problem
        and not has_hierarchical_blocks
        and is_complex_clean2d_graph(graph, selected)
    ):
        return _run_local_graph_clean2d_engine(graph, selected, before, mode, target)

    candidates = generate_clean2d_candidates(
        graph,
        selected,
        mode=mode,
        target_bond_length=target_bond_length,
        avoid_hashes=avoid_hashes,
        seed=seed,
        rdkit_timeout_s=rdkit_timeout_s,
    )
    return rank_clean2d_candidates(
        graph,
        candidates,
        before,
        selected,
        mode=mode,
        target_bond_length=target_bond_length,
        avoid_hashes=avoid_hashes,
    )


def _run_local_graph_clean2d_engine(
    graph: MolGraph,
    selected: set[int],
    before: dict[int, tuple[float, float]],
    mode: Clean2DMode,
    target: float,
) -> Clean2DResult:
    local_mode = (
        LocalClean2DMode.PUBLICATION
        if mode == Clean2DMode.PUBLICATION
        else LocalClean2DMode.QUICK
    )
    local_result = local_graph_clean2d(
        graph,
        selected,
        target_bond_length=target,
        mode=local_mode,
    )
    bonds = _selected_structural_bonds(graph, selected)
    metadata = {
        "local_graph_clean2d": True,
        **local_result.report.as_dict(),
    }

    if local_result.ok:
        after = _complete_coords(local_result.coords, before, selected)
        report = evaluate_clean2d_layout(
            selected,
            bonds,
            before,
            after,
            target,
            is_cyclic=has_cycles(selected, bonds),
        )
        quality = classify_clean2d_layout_quality(
            graph,
            selected,
            coords=after,
            target_bond_length=target,
        )
        metadata.update(
            {
                "quality_class": quality.quality_class,
                "quality_reason": quality.reason,
                "length_rms_error": quality.length_rms_error,
                "length_max_error": quality.length_max_error,
                "angle_rms_deviation": quality.angle_rms_deviation,
                "angle_max_deviation": quality.angle_max_deviation,
                "severe_angle_count": quality.severe_angle_count,
                "crossings": quality.crossings,
                "min_nonbonded_distance": quality.min_nonbonded_distance,
                "min_ring_degeneracy": quality.min_ring_degeneracy,
                "exocyclic_ring_angle_min": quality.exocyclic_ring_angle_min,
                "bad_exocyclic_ring_count": quality.bad_exocyclic_ring_count,
                "visual_score": quality.visual_score,
            }
        )
        candidate = Clean2DCandidate(
            source="local_graph",
            coords=after,
            message=local_result.report.message or "Estructura 2D limpiada localmente",
            score=_visual_quality_score(graph, selected, bonds, before, after, target, mode),
            novelty=_mean_displacement(before, after, selected),
            report=report,
            geometry_hash=clean2d_geometry_hash(graph, after, selected),
            metadata=metadata,
        )
        return Clean2DResult(
            mode=mode,
            atom_ids=selected,
            before_coords=before,
            candidates=(candidate,),
            selected=candidate,
            message=candidate.message,
        )

    rejected = Clean2DCandidate(
        source="local_graph",
        coords=before,
        message=local_result.report.message,
        rejected=True,
        rejection_reason=local_result.report.reason or "no_safe_local_repair",
        geometry_hash=clean2d_geometry_hash(graph, before, selected) if before else "",
        metadata=metadata,
    )
    return Clean2DResult(
        mode=mode,
        atom_ids=selected,
        before_coords=before,
        candidates=(),
        selected=None,
        rejected=(rejected,),
        message=local_result.report.message,
        reason=local_result.report.reason or "no_safe_local_repair",
    )


def generate_clean2d_candidates(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    mode: Clean2DMode | str = Clean2DMode.QUICK,
    target_bond_length: float = 42.0,
    avoid_hashes: set[str] | None = None,
    seed: int | None = None,
    rdkit_timeout_s: float = 8.0,
) -> list[Clean2DCandidate]:
    mode = _coerce_mode(mode)
    selected = _normalize_atom_ids(graph, atom_ids)
    if not selected:
        return []

    before = _coords_for_atoms(graph, selected)
    bonds = _selected_structural_bonds(graph, selected)
    target = max(8.0, float(target_bond_length or 42.0))
    layer_model = build_multilayer_chemical_graph(graph, selected)
    baseline_constraint_error = _interaction_constraint_error(layer_model.layout_constraint_graph.constraints, before)
    block_layout_candidate = _candidate_from_block_layout_signals(graph, selected, before, bonds, target)
    block_candidate = _candidate_from_block_constraints(graph, selected, before, bonds, target)
    motif_candidate = _candidate_from_motif_constraints(graph, selected, before, bonds, target)
    candidates: list[Clean2DCandidate] = []
    if block_layout_candidate is not None:
        candidates.append(block_layout_candidate)

    current = _normalized_candidate_coords(before, before, bonds, target)
    candidates.append(
        Clean2DCandidate(
            source="current",
            coords=current,
            message="Estructura 2D ya estaba limpia",
            metadata={
                "interaction_constraint_error": baseline_constraint_error,
                "interaction_constraint_count": _interaction_constraint_count(
                    layer_model.layout_constraint_graph.constraints
                ),
            },
        )
    )

    if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}:
        candidates.extend(
            candidate
            for candidate in (
                block_candidate,
                motif_candidate,
                _candidate_from_rdkit_isolated(graph, selected, before, bonds, target, rdkit_timeout_s),
                _candidate_from_rdkit_direct(graph, selected, before, bonds, target),
                _candidate_from_internal_templates(graph, selected, before, bonds, target),
                _candidate_from_v2(graph, selected, before, bonds, target, mode),
                _candidate_from_structure_preserving_fallback(selected, before, bonds, target),
            )
            if candidate is not None
        )
    else:
        candidates.extend(_propose_candidates(graph, selected, before, bonds, target, seed=seed))
        candidates.extend(
            candidate
            for candidate in (
                block_layout_candidate,
                block_candidate,
                motif_candidate,
                _candidate_from_rdkit_isolated(graph, selected, before, bonds, target, rdkit_timeout_s),
                _candidate_from_internal_templates(graph, selected, before, bonds, target),
            )
            if candidate is not None
        )

    return _deduplicate_candidates(graph, selected, candidates, avoid_hashes=avoid_hashes)


def rank_clean2d_candidates(
    graph: MolGraph,
    candidates: Iterable[Clean2DCandidate],
    before_coords: dict[int, tuple[float, float]],
    atom_ids: Iterable[int],
    *,
    mode: Clean2DMode | str = Clean2DMode.QUICK,
    target_bond_length: float = 42.0,
    avoid_hashes: set[str] | None = None,
) -> Clean2DResult:
    mode = _coerce_mode(mode)
    selected = set(atom_ids)
    bonds = _selected_structural_bonds(graph, selected)
    target = max(8.0, float(target_bond_length or 42.0))
    cyclic = has_cycles(selected, bonds)
    baseline_score = _visual_quality_score(graph, selected, bonds, before_coords, before_coords, target, mode)
    quality = classify_clean2d_layout_quality(
        graph,
        selected,
        coords=before_coords,
        target_bond_length=target,
    )
    baseline_needs_work = quality.quality_class != "good"
    baseline_bad = quality.quality_class == "needs_rebuild"
    avoid_hashes = set(avoid_hashes or set())

    if mode == Clean2DMode.PROPOSE and quality.quality_class != "good":
        rejected = []
        reason = "geometria_base_requiere_optimizacion"
        for candidate in candidates:
            after = _complete_coords(candidate.coords, before_coords, selected) if candidate.coords else {}
            rejected.append(
                _replace_candidate(
                    candidate,
                    coords=after,
                    rejected=True,
                    rejection_reason=reason,
                    geometry_hash=(
                        clean2d_geometry_hash(graph, after, selected)
                        if after
                        else candidate.geometry_hash
                    ),
                    metadata={**candidate.metadata, "quality_class": quality.quality_class},
                )
            )
        return Clean2DResult(
            mode=mode,
            atom_ids=selected,
            before_coords=before_coords,
            candidates=(),
            selected=None,
            rejected=tuple(rejected),
            message="La estructura debe optimizarse antes de proponer conformeros 2D",
            reason=reason,
        )

    accepted: list[Clean2DCandidate] = []
    rejected: list[Clean2DCandidate] = []
    for candidate in candidates:
        if candidate.rejected:
            rejected.append(candidate)
            continue
        after = _complete_coords(candidate.coords, before_coords, selected)
        geometry_hash = clean2d_geometry_hash(graph, after, selected)
        try:
            assert_clean2d_invariants(
                graph,
                graph,
                before_coords,
                after,
                atom_ids=selected,
            )
        except Clean2DInvariantError as exc:
            rejected.append(
                _replace_candidate(
                    candidate,
                    coords=after,
                    rejected=True,
                    rejection_reason=str(exc),
                    geometry_hash=geometry_hash,
                )
            )
            continue

        report = evaluate_clean2d_layout(
            selected,
            bonds,
            before_coords,
            after,
            target,
            is_cyclic=cyclic,
        )
        candidate_quality = classify_clean2d_layout_quality(
            graph,
            selected,
            coords=after,
            target_bond_length=target,
        )
        candidate_metadata = {
            **candidate.metadata,
            "quality_class": candidate_quality.quality_class,
            "quality_reason": candidate_quality.reason,
            "length_rms_error": candidate_quality.length_rms_error,
            "length_max_error": candidate_quality.length_max_error,
            "angle_rms_deviation": candidate_quality.angle_rms_deviation,
            "angle_max_deviation": candidate_quality.angle_max_deviation,
            "severe_angle_count": candidate_quality.severe_angle_count,
            "crossings": candidate_quality.crossings,
            "min_nonbonded_distance": candidate_quality.min_nonbonded_distance,
            "min_ring_degeneracy": candidate_quality.min_ring_degeneracy,
            "exocyclic_ring_angle_min": candidate_quality.exocyclic_ring_angle_min,
            "bad_exocyclic_ring_count": candidate_quality.bad_exocyclic_ring_count,
            "visual_score": candidate_quality.visual_score,
        }
        if candidate_metadata.get("interaction_constraint_count"):
            candidate_metadata["interaction_constraint_error"] = _interaction_constraint_error(
                build_multilayer_chemical_graph(graph, selected).layout_constraint_graph.constraints,
                after,
            )
        safety_mode = _safety_mode_for_candidate(mode, candidate, baseline_bad)
        safe = is_clean2d_candidate_safe(report, safety_mode)
        if candidate.source == "block_layout" and baseline_bad and not candidate.metadata.get("block_layout_rejected"):
            safe = True
        if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION} and candidate.source == "current":
            if quality.quality_class == "needs_rebuild":
                safe = False
                report.rejection_reason = "geometria_actual_necesita_reconstruccion"
            elif quality.quality_class == "needs_polish":
                safe = False
                report.rejection_reason = "geometria_actual_requiere_optimizacion"
        if (
            mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}
            and candidate.source != "block_layout"
            and candidate_quality.quality_class == "needs_rebuild"
        ):
            safe = False
            report.rejection_reason = "candidato_no_canonicaliza_geometria"
        if (
            mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}
            and baseline_needs_work
            and candidate_quality.quality_class == "needs_polish"
            and not _candidate_substantially_improves_layout(quality, candidate_quality)
        ):
            safe = False
            report.rejection_reason = "candidato_no_mejora_suficiente"
        if mode == Clean2DMode.PROPOSE and candidate.source == "current":
            safe = False
            report.rejection_reason = "propose_requiere_geometria_alternativa"
        if (
            mode == Clean2DMode.PROPOSE
            and not baseline_bad
            and geometry_hash in avoid_hashes
            and candidate.source != "current"
        ):
            safe = False
            report.rejection_reason = "geometria_repetida"
        if not safe:
            rejected.append(
                _replace_candidate(
                    candidate,
                    coords=after,
                    rejected=True,
                    rejection_reason=report.rejection_reason or "candidato_rechazado",
                    report=report,
                    geometry_hash=geometry_hash,
                    metadata=candidate_metadata,
                )
            )
            continue

        novelty = _mean_displacement(before_coords, after, selected)
        score = _visual_quality_score(graph, selected, bonds, before_coords, after, target, mode)
        if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION} and baseline_needs_work:
            score += _quality_rank(candidate_quality.quality_class) * target * 10000.0
            score += candidate_quality.visual_score
            score += candidate_quality.length_rms_error * target * 120.0
            score += candidate_quality.angle_rms_deviation * 12.0
            score += candidate_quality.crossings * target * 1000.0
            score += candidate_quality.bad_exocyclic_ring_count * target * 500.0
            if candidate_quality.exocyclic_ring_angle_min < math.inf:
                score += max(0.0, 125.0 - candidate_quality.exocyclic_ring_angle_min) * target * 2.0
            if candidate_quality.min_nonbonded_distance < math.inf:
                score -= min(candidate_quality.min_nonbonded_distance, target * 2.0) * 0.05
            score += novelty * 0.02
        elif mode == Clean2DMode.QUICK and candidate.source != "current":
            score += novelty * 12.0
        elif mode == Clean2DMode.QUICK:
            score += novelty * 0.20
        elif mode == Clean2DMode.PUBLICATION:
            score += novelty * 0.05
        else:
            novelty_bonus = min(target * 2.0, novelty)
            score -= novelty_bonus * 3.0
            if novelty < target * 0.05:
                score += target * 80.0

        message = candidate.message
        if (
            mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}
            and baseline_needs_work
            and candidate_quality.quality_class == "needs_polish"
        ):
            message = "Estructura 2D parcialmente optimizada; requiere otra pasada"
        elif (
            mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}
            and baseline_needs_work
            and candidate_quality.quality_class == "good"
            and candidate.source != "current"
        ):
            message = "Estructura 2D optimizada"

        accepted.append(
            _replace_candidate(
                candidate,
                coords=after,
                message=message,
                score=score,
                novelty=novelty,
                report=report,
                geometry_hash=geometry_hash,
                metadata=candidate_metadata,
            )
        )

    if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION} and baseline_needs_work:
        accepted.sort(key=lambda item: _canonicalization_sort_key(item))
    elif mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}:
        accepted = _rank_good_baseline_candidates_for_canonicalization(accepted, target)
    elif mode == Clean2DMode.PROPOSE:
        accepted.sort(
            key=lambda item: (
                round(item.score, 9),
                _propose_source_priority(item.source),
                item.score,
                item.source,
            )
        )
    else:
        accepted.sort(key=lambda item: (item.score, _source_priority(item.source), item.source))
    selected_candidate = accepted[0] if accepted else None
    message = selected_candidate.message if selected_candidate is not None else ""
    reason = ""
    if selected_candidate is None:
        reason = rejected[0].rejection_reason if rejected else "sin_candidatos"
        message = f"Limpieza 2D omitida: {reason}"

    return Clean2DResult(
        mode=mode,
        atom_ids=selected,
        before_coords=before_coords,
        candidates=tuple(accepted),
        selected=selected_candidate,
        rejected=tuple(rejected),
        message=message,
        reason=reason,
)


def classify_clean2d_layout_quality(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    coords: dict[int, tuple[float, float]] | None = None,
    target_bond_length: float = 42.0,
) -> Clean2DLayoutQualityReport:
    selected = _normalize_atom_ids(graph, atom_ids)
    target = max(8.0, float(target_bond_length or 42.0))
    before = _coords_for_atoms(graph, selected) if coords is None else _complete_coords(coords, _coords_for_atoms(graph, selected), selected)
    if not selected:
        return Clean2DLayoutQualityReport("good", reason="sin_atomos")
    missing = selected - set(before)
    if missing:
        return Clean2DLayoutQualityReport("needs_rebuild", reason="faltan_coordenadas")
    if any(not (math.isfinite(x) and math.isfinite(y)) for x, y in before.values()):
        return Clean2DLayoutQualityReport("needs_rebuild", reason="coordenadas_no_finitas")

    bonds = _selected_structural_bonds(graph, selected)
    if not bonds:
        return Clean2DLayoutQualityReport("good", reason="sin_enlaces")

    length_errors = []
    for bond in bonds:
        if bond.a1_id not in before or bond.a2_id not in before:
            return Clean2DLayoutQualityReport("needs_rebuild", reason="faltan_coordenadas_enlace")
        observed = _distance(before[bond.a1_id], before[bond.a2_id])
        desired = _target_length_for_bond(bond, target)
        length_errors.append(abs(observed - desired) / max(target, 1e-6))
    length_rms = math.sqrt(sum(err * err for err in length_errors) / len(length_errors))
    length_max = max(length_errors)

    angle_stats = _strict_angle_deviation_stats(graph, selected, bonds, before)
    crossings = _count_crossings(before, bonds)
    nonbonded = min_nonbonded_distance(before, bonds, selected)
    ring_scores = [
        ring_degeneracy_score(before, set(ring))
        for ring in _cycle_basis_ordered(selected, bonds, max_size=8)
    ]
    min_ring = min(ring_scores) if ring_scores else math.inf
    exocyclic_stats = _exocyclic_ring_orientation_stats(selected, bonds, before)
    visual_score = _visual_quality_score(graph, selected, bonds, before, before, target, Clean2DMode.QUICK)

    base_report = Clean2DLayoutQualityReport(
        "good",
        length_rms_error=length_rms,
        length_max_error=length_max,
        angle_rms_deviation=angle_stats["rms"],
        angle_max_deviation=angle_stats["max"],
        angle_penalty=angle_stats["penalty"],
        severe_angle_count=int(angle_stats["severe_count"]),
        crossings=crossings,
        min_nonbonded_distance=nonbonded,
        min_ring_degeneracy=min_ring,
        exocyclic_ring_angle_min=exocyclic_stats["min_center_angle"],
        bad_exocyclic_ring_count=int(exocyclic_stats["bad_count"]),
        visual_score=visual_score,
    )

    if _layout_needs_rebuild(selected, bonds, before, target):
        return _replace_quality(base_report, quality_class="needs_rebuild", reason="fallos_geometricos_graves")
    if length_max > 0.42 or length_rms > 0.22:
        return _replace_quality(base_report, quality_class="needs_rebuild", reason="longitudes_fuera_de_rango")
    if crossings > 0:
        return _replace_quality(base_report, quality_class="needs_rebuild", reason="cruces_enlaces")
    if nonbonded < target * 0.35:
        return _replace_quality(base_report, quality_class="needs_rebuild", reason="colisiones_no_enlazadas")
    if min_ring != math.inf and min_ring < 0.12:
        return _replace_quality(base_report, quality_class="needs_rebuild", reason="anillo_degenerado")
    if exocyclic_stats["bad_count"] > 0 and (crossings > 0 or nonbonded < target * 0.35):
        return _replace_quality(
            base_report,
            quality_class="needs_rebuild",
            reason="orientacion_anillo_exociclico_con_colision",
        )

    if length_max > 0.14 or length_rms > 0.055:
        return _replace_quality(base_report, quality_class="needs_polish", reason="longitudes_suboptimas")
    if angle_stats["max"] > 34.0 or angle_stats["rms"] > 18.0 or angle_stats["severe_count"] > 0:
        return _replace_quality(base_report, quality_class="needs_polish", reason="angulos_suboptimos")
    if nonbonded < target * 0.50:
        return _replace_quality(base_report, quality_class="needs_polish", reason="distancias_no_enlazadas_estrechas")
    if min_ring != math.inf and min_ring < 0.25:
        return _replace_quality(base_report, quality_class="needs_polish", reason="anillo_suboptimo")
    if exocyclic_stats["bad_count"] > 0:
        return _replace_quality(
            base_report,
            quality_class="needs_polish",
            reason="orientacion_anillo_exociclico_suboptima",
        )
    return base_report


def has_intramolecular_block_layout_problem(
    block_graph: object,
    quality_report: Clean2DLayoutQualityReport | None,
) -> bool:
    if has_hierarchical_block_layout_signals(block_graph):
        return True
    blocks = tuple(getattr(block_graph, "blocks", ()) or ())
    if not blocks:
        return False
    counts: dict[BlockKind, int] = {}
    for block in blocks:
        kind = getattr(block, "kind", None)
        if isinstance(kind, BlockKind):
            counts[kind] = counts.get(kind, 0) + 1
    structural_trigger = any(
        counts.get(kind, 0) > 0
        for kind in (
            BlockKind.MACROCYCLE,
            BlockKind.CYCLOPHANE,
            BlockKind.FUSED_SYSTEM,
            BlockKind.INTRAMOLECULAR_BRIDGE,
            BlockKind.INTERNAL_CAVITY,
        )
    )
    if structural_trigger:
        return True
    if counts.get(BlockKind.AROMATIC_RING, 0) >= 3:
        return True
    if counts.get(BlockKind.STEREO_CENTER, 0) > 0 and counts.get(BlockKind.AROMATIC_RING, 0) > 0:
        return True
    if quality_report is None:
        return False
    return quality_report.quality_class != "good" and (
        counts.get(BlockKind.AROMATIC_RING, 0) >= 2
        or counts.get(BlockKind.LINKER, 0) > 0
        or counts.get(BlockKind.TERMINAL_SUBSTITUENT, 0) > 0
    )


def has_hierarchical_block_layout_signals(block_graph: object) -> bool:
    blocks = tuple(getattr(block_graph, "blocks", ()) or ())
    if not blocks:
        return False
    counts: dict[BlockKind, int] = {}
    block_by_id = {}
    for block in blocks:
        block_id = int(getattr(block, "id", 0) or 0)
        if block_id:
            block_by_id[block_id] = block
        kind = getattr(block, "kind", None)
        if isinstance(kind, BlockKind):
            counts[kind] = counts.get(kind, 0) + 1
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
    if counts.get(BlockKind.STEREO_CENTER, 0) >= 2:
        return True

    rigid = {
        BlockKind.AROMATIC_RING,
        BlockKind.FUSED_SYSTEM,
        BlockKind.MACROCYCLE,
        BlockKind.CYCLOPHANE,
    }
    linker_edges = 0
    for edge in getattr(block_graph, "edges", ()) or ():
        if getattr(edge, "kind", None) != BlockEdgeKind.LINKER:
            continue
        block_ids = tuple(getattr(edge, "block_ids", ()) or ())
        if len(block_ids) != 2:
            continue
        left = block_by_id.get(int(block_ids[0]))
        right = block_by_id.get(int(block_ids[1]))
        if left is None or right is None:
            continue
        if getattr(left, "kind", None) in rigid or getattr(right, "kind", None) in rigid:
            linker_edges += 1
    return linker_edges >= 2


def _replace_quality(report: Clean2DLayoutQualityReport, **updates: Any) -> Clean2DLayoutQualityReport:
    data = {
        "quality_class": report.quality_class,
        "reason": report.reason,
        "length_rms_error": report.length_rms_error,
        "length_max_error": report.length_max_error,
        "angle_rms_deviation": report.angle_rms_deviation,
        "angle_max_deviation": report.angle_max_deviation,
        "angle_penalty": report.angle_penalty,
        "severe_angle_count": report.severe_angle_count,
        "crossings": report.crossings,
        "min_nonbonded_distance": report.min_nonbonded_distance,
        "min_ring_degeneracy": report.min_ring_degeneracy,
        "exocyclic_ring_angle_min": report.exocyclic_ring_angle_min,
        "bad_exocyclic_ring_count": report.bad_exocyclic_ring_count,
        "visual_score": report.visual_score,
    }
    data.update(updates)
    return Clean2DLayoutQualityReport(**data)


def summarize_clean2d_candidates(result: Clean2DResult) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate in (*result.candidates, *result.rejected):
        metadata = candidate.metadata
        rows.append(
            {
                "source": candidate.source,
                "rejected": candidate.rejected,
                "reason": candidate.rejection_reason,
                "score": candidate.score,
                "novelty": candidate.novelty,
                "quality_class": metadata.get("quality_class", ""),
                "quality_reason": metadata.get("quality_reason", ""),
                "length_rms_error": metadata.get("length_rms_error", 0.0),
                "length_max_error": metadata.get("length_max_error", 0.0),
                "angle_rms_deviation": metadata.get("angle_rms_deviation", 0.0),
                "angle_max_deviation": metadata.get("angle_max_deviation", 0.0),
                "severe_angle_count": metadata.get("severe_angle_count", 0),
                "crossings": metadata.get("crossings", 0),
                "min_nonbonded_distance": metadata.get("min_nonbonded_distance", math.inf),
                "min_ring_degeneracy": metadata.get("min_ring_degeneracy", math.inf),
                "exocyclic_ring_angle_min": metadata.get("exocyclic_ring_angle_min", math.inf),
                "bad_exocyclic_ring_count": metadata.get("bad_exocyclic_ring_count", 0),
                "visual_score": metadata.get("visual_score", 0.0),
                "geometry_hash": candidate.geometry_hash,
            }
        )
    return rows


def _quality_rank(quality_class: str) -> int:
    return {
        "good": 0,
        "needs_polish": 1,
        "needs_rebuild": 2,
    }.get(str(quality_class or ""), 9)


def _candidate_substantially_improves_layout(
    before: Clean2DLayoutQualityReport,
    after: Clean2DLayoutQualityReport,
) -> bool:
    if after.quality_class == "good":
        return True
    if after.quality_class == "needs_rebuild":
        return False
    no_new_crossings = after.crossings <= before.crossings
    no_new_collision = (
        after.min_nonbonded_distance >= before.min_nonbonded_distance - 1e-6
        or before.min_nonbonded_distance == math.inf
    )
    no_ring_worse = (
        after.min_ring_degeneracy >= before.min_ring_degeneracy - 1e-6
        or before.min_ring_degeneracy == math.inf
    )
    if not (no_new_crossings and no_new_collision and no_ring_worse):
        return False
    if before.visual_score > 1e-6 and after.visual_score < before.visual_score * 0.75:
        return True
    length_better = (
        after.length_rms_error < before.length_rms_error * 0.85
        or after.length_rms_error + 0.01 < before.length_rms_error
    )
    angle_better = (
        after.angle_rms_deviation < before.angle_rms_deviation * 0.85
        or after.angle_rms_deviation + 2.0 < before.angle_rms_deviation
    )
    return length_better and angle_better


def _rank_good_baseline_candidates_for_canonicalization(
    accepted: list[Clean2DCandidate],
    target: float,
) -> list[Clean2DCandidate]:
    if not accepted:
        return []
    current = next((candidate for candidate in accepted if candidate.source == "current"), None)
    if current is None:
        return sorted(accepted, key=lambda item: (item.score, _source_priority(item.source), item.source))
    if not _current_has_material_canonical_defect(current.metadata, target):
        current = _replace_candidate(
            current,
            metadata={**current.metadata, "current_canonical_enough": True},
        )
        remainder = [candidate for candidate in accepted if candidate.source != "current"]
        return [current] + sorted(remainder, key=lambda item: (item.score, _source_priority(item.source), item.source))

    improvements = [
        _with_canonical_delta(current, candidate, target)
        for candidate in accepted
        if candidate.source != "current"
    ]
    improvements = [
        candidate
        for candidate in improvements
        if bool(candidate.metadata.get("canonical_improvement_over_current", False))
    ]
    if improvements:
        improvements.sort(key=lambda item: _canonicalization_sort_key(item))
        chosen = _replace_candidate(improvements[0], message="Estructura 2D optimizada")
        remainder = [
            candidate
            for candidate in accepted
            if candidate.geometry_hash != chosen.geometry_hash or candidate.source != chosen.source
        ]
        return [chosen] + sorted(remainder, key=lambda item: (item.score, _source_priority(item.source), item.source))

    current = _replace_candidate(
        current,
        metadata={**current.metadata, "current_canonical_enough": True},
    )
    remainder = [candidate for candidate in accepted if candidate.source != "current"]
    return [current] + sorted(remainder, key=lambda item: (item.score, _source_priority(item.source), item.source))


def _current_has_material_canonical_defect(metadata: dict[str, Any], target: float) -> bool:
    visual_score = float(metadata.get("visual_score", 0.0) or 0.0)
    if visual_score > target * 20.0:
        return True
    if float(metadata.get("bad_exocyclic_ring_count", 0) or 0) > 0:
        return True
    if float(metadata.get("crossings", 0) or 0) > 0:
        return True
    if float(metadata.get("severe_angle_count", 0) or 0) > 0:
        return True
    if float(metadata.get("interaction_constraint_error", 0.0) or 0.0) > target * 0.25:
        return True
    min_nonbonded = float(metadata.get("min_nonbonded_distance", math.inf) or math.inf)
    if min_nonbonded < target * 0.50:
        return True
    min_ring = float(metadata.get("min_ring_degeneracy", math.inf) or math.inf)
    return min_ring < 0.25


def _with_canonical_delta(
    current: Clean2DCandidate,
    candidate: Clean2DCandidate,
    target: float,
) -> Clean2DCandidate:
    delta = canonicalization_delta(current, candidate, target)
    return _replace_candidate(candidate, metadata={**candidate.metadata, **delta})


def canonicalization_delta(
    current: Clean2DCandidate,
    candidate: Clean2DCandidate,
    target: float,
) -> dict[str, Any]:
    current_meta = current.metadata
    candidate_meta = candidate.metadata
    novelty = float(candidate.novelty or 0.0)
    current_visual = float(current_meta.get("visual_score", current.score) or current.score or 0.0)
    candidate_visual = float(candidate_meta.get("visual_score", candidate.score) or candidate.score or math.inf)
    current_quality = str(current_meta.get("quality_class", ""))
    candidate_quality = str(candidate_meta.get("quality_class", ""))
    geometry_equivalent = (
        bool(current.geometry_hash)
        and bool(candidate.geometry_hash)
        and current.geometry_hash == candidate.geometry_hash
    ) or novelty < target * 0.05
    metric_better = _candidate_improves_important_metric(current_meta, candidate_meta, target)
    visual_margin = max(0.5, target * 0.02)
    visual_better = candidate_visual + visual_margin < current_visual
    visual_close_for_dirty_layout = (
        current_visual > target * 20.0
        and candidate_visual <= max(current_visual, 1e-6) * 1.05
    )
    canonical_improvement = (
        current_quality == "good"
        and candidate_quality == "good"
        and not geometry_equivalent
        and (visual_better or visual_close_for_dirty_layout or metric_better)
    )
    return {
        "current_canonical_enough": False,
        "canonical_improvement_over_current": canonical_improvement,
        "canonical_geometry_equivalent": geometry_equivalent,
        "canonical_visual_delta": current_visual - candidate_visual,
        "canonical_metric_better": metric_better,
    }


def _candidate_improves_important_metric(
    current_meta: dict[str, Any],
    candidate_meta: dict[str, Any],
    target: float,
) -> bool:
    if float(candidate_meta.get("crossings", 0) or 0) < float(current_meta.get("crossings", 0) or 0):
        return True
    if float(candidate_meta.get("bad_exocyclic_ring_count", 0) or 0) < float(
        current_meta.get("bad_exocyclic_ring_count", 0) or 0
    ):
        return True
    if float(candidate_meta.get("severe_angle_count", 0) or 0) < float(current_meta.get("severe_angle_count", 0) or 0):
        return True

    current_angle = float(current_meta.get("angle_rms_deviation", 0.0) or 0.0)
    candidate_angle = float(candidate_meta.get("angle_rms_deviation", 0.0) or 0.0)
    if current_angle > 3.0 and (candidate_angle < current_angle * 0.85 or candidate_angle + 2.0 < current_angle):
        return True

    current_length = float(current_meta.get("length_rms_error", 0.0) or 0.0)
    candidate_length = float(candidate_meta.get("length_rms_error", 0.0) or 0.0)
    if current_length > 0.055 and (candidate_length < current_length * 0.85 or candidate_length + 0.01 < current_length):
        return True

    current_ring = float(current_meta.get("min_ring_degeneracy", math.inf) or math.inf)
    candidate_ring = float(candidate_meta.get("min_ring_degeneracy", math.inf) or math.inf)
    if current_ring < 0.35 and candidate_ring > current_ring + 0.05:
        return True

    current_exocyclic = float(current_meta.get("exocyclic_ring_angle_min", math.inf) or math.inf)
    candidate_exocyclic = float(candidate_meta.get("exocyclic_ring_angle_min", math.inf) or math.inf)
    if current_exocyclic < 125.0 and candidate_exocyclic > current_exocyclic + 15.0:
        return True

    current_interaction = float(current_meta.get("interaction_constraint_error", 0.0) or 0.0)
    candidate_interaction = float(candidate_meta.get("interaction_constraint_error", 0.0) or 0.0)
    if current_interaction > target * 0.20 and candidate_interaction < current_interaction * 0.55:
        return True

    current_nonbonded = float(current_meta.get("min_nonbonded_distance", math.inf) or math.inf)
    candidate_nonbonded = float(candidate_meta.get("min_nonbonded_distance", math.inf) or math.inf)
    return current_nonbonded < target * 0.50 and candidate_nonbonded > current_nonbonded + target * 0.10


def _canonicalization_sort_key(candidate: Clean2DCandidate) -> tuple[float, ...]:
    metadata = candidate.metadata
    quality_class = str(metadata.get("quality_class", ""))
    min_nonbonded = float(metadata.get("min_nonbonded_distance", math.inf))
    min_ring = float(metadata.get("min_ring_degeneracy", math.inf))
    return (
        float(_quality_rank(quality_class)),
        float(metadata.get("visual_score", candidate.score) or candidate.score),
        float(metadata.get("length_rms_error", 0.0) or 0.0),
        float(metadata.get("angle_rms_deviation", 0.0) or 0.0),
        float(metadata.get("crossings", 0) or 0),
        float(metadata.get("length_max_error", 0.0) or 0.0),
        float(metadata.get("angle_max_deviation", 0.0) or 0.0),
        float(metadata.get("severe_angle_count", 0) or 0),
        float(metadata.get("bad_exocyclic_ring_count", 0) or 0),
        -float(metadata.get("exocyclic_ring_angle_min", math.inf) or math.inf),
        -min(min_nonbonded, 1_000_000.0),
        -min(min_ring, 1_000_000.0),
        float(_source_priority(candidate.source)),
    )


def capture_clean2d_snapshot(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
) -> Clean2DGraphSnapshot:
    selected = _normalize_atom_ids(graph, atom_ids)
    bonds = [
        bond
        for bond in graph.bonds.values()
        if bond.a1_id in selected and bond.a2_id in selected
    ]
    atom_data = []
    for atom_id in sorted(selected):
        atom = graph.atoms[atom_id]
        atom_data.append(
            (
                atom_id,
                (
                    atom.element,
                    int(getattr(atom, "charge", 0) or 0),
                    getattr(atom, "isotope", None),
                    int(getattr(atom, "radical_electrons", 0) or 0),
                    getattr(atom, "stereo_cip", None),
                    getattr(atom, "stereo_axial", None),
                    getattr(atom, "stereo_helical", None),
                    getattr(atom, "stereo_si_re", None),
                    getattr(atom, "explicit_h", None),
                    getattr(atom, "group_h_cap", None),
                    getattr(atom, "mapping", None),
                    bool(getattr(atom, "is_query", False)),
                    bool(getattr(atom, "is_explicit", False)),
                    bool(getattr(atom, "no_implicit", False)),
                    tuple(getattr(atom, "r_group_substituents", ()) or ()),
                ),
            )
        )

    bond_data = []
    for bond in sorted(bonds, key=lambda item: item.id):
        bond_data.append(
            (
                int(bond.id),
                (
                    int(bond.a1_id),
                    int(bond.a2_id),
                    int(getattr(bond, "order", 1) or 1),
                    getattr(getattr(bond, "style", None), "value", getattr(bond, "style", None)),
                    getattr(getattr(bond, "stereo", None), "value", getattr(bond, "stereo", None)),
                    getattr(bond, "stereo_ez", None),
                    getattr(bond, "stereo_axial", None),
                    getattr(bond, "stereo_endo_exo", None),
                    getattr(bond, "stereo_helical", None),
                    bool(getattr(bond, "is_aromatic", False)),
                    getattr(bond, "display_order", None),
                    bool(getattr(bond, "is_query", False)),
                    getattr(bond, "ring_id", None),
                    getattr(bond, "donor_atom_id", None),
                    getattr(bond, "pi_offset_sign", None),
                    getattr(bond, "interaction_kind", None),
                ),
            )
        )

    return Clean2DGraphSnapshot(
        atom_ids=frozenset(selected),
        bond_ids=frozenset(int(bond.id) for bond in bonds),
        atom_data=tuple(atom_data),
        bond_data=tuple(bond_data),
        components=_components_for_snapshot(selected, bonds),
    )


def assert_clean2d_invariants(
    before_graph: MolGraph | Clean2DGraphSnapshot,
    after_graph: MolGraph | Clean2DGraphSnapshot,
    before_coords: dict[int, tuple[float, float]],
    after_coords: dict[int, tuple[float, float]],
    *,
    atom_ids: Iterable[int] | None = None,
) -> None:
    before_snapshot = (
        before_graph
        if isinstance(before_graph, Clean2DGraphSnapshot)
        else capture_clean2d_snapshot(before_graph, atom_ids)
    )
    after_snapshot = (
        after_graph
        if isinstance(after_graph, Clean2DGraphSnapshot)
        else capture_clean2d_snapshot(after_graph, atom_ids)
    )
    if before_snapshot.atom_ids != after_snapshot.atom_ids:
        raise Clean2DInvariantError("cambio_conjunto_atomos")
    if before_snapshot.bond_ids != after_snapshot.bond_ids:
        raise Clean2DInvariantError("cambio_conjunto_enlaces")
    if before_snapshot.atom_data != after_snapshot.atom_data:
        raise Clean2DInvariantError("cambio_metadatos_atomos")
    if before_snapshot.bond_data != after_snapshot.bond_data:
        raise Clean2DInvariantError("cambio_metadatos_enlaces")
    if before_snapshot.components != after_snapshot.components:
        raise Clean2DInvariantError("cambio_componentes_conectadas")

    selected = set(atom_ids) if atom_ids is not None else set(before_snapshot.atom_ids)
    if selected != set(before_snapshot.atom_ids):
        selected &= set(before_snapshot.atom_ids)
    for atom_id in selected:
        if atom_id not in before_coords:
            raise Clean2DInvariantError(f"coordenada_inicial_faltante:{atom_id}")
        if atom_id not in after_coords:
            raise Clean2DInvariantError(f"coordenada_final_faltante:{atom_id}")
        x, y = after_coords[atom_id]
        if not (math.isfinite(float(x)) and math.isfinite(float(y))):
            raise Clean2DInvariantError(f"coordenada_no_finita:{atom_id}")


def clean2d_geometry_hash(
    graph: MolGraph,
    coords: dict[int, tuple[float, float]],
    atom_ids: Iterable[int] | None = None,
    *,
    precision: float = 0.25,
) -> str:
    selected = _normalize_atom_ids(graph, atom_ids)
    bonds = _selected_structural_bonds(graph, selected)
    center = _center({atom_id: coords[atom_id] for atom_id in selected if atom_id in coords})
    avg_len = _average_bond_length(coords, bonds) or 1.0
    scale = max(avg_len, 1e-6)
    rounded_coords = []
    for atom_id in sorted(selected):
        x, y = coords.get(atom_id, (0.0, 0.0))
        nx = (float(x) - center[0]) / scale
        ny = (float(y) - center[1]) / scale
        rounded_coords.append((atom_id, round(nx / precision), round(ny / precision)))
    payload = repr(
        (
            tuple(
                (
                    atom_id,
                    graph.atoms[atom_id].element,
                    int(getattr(graph.atoms[atom_id], "charge", 0) or 0),
                )
                for atom_id in sorted(selected)
            ),
            tuple(
                (
                    min(bond.a1_id, bond.a2_id),
                    max(bond.a1_id, bond.a2_id),
                    int(getattr(bond, "order", 1) or 1),
                    bool(getattr(bond, "is_aromatic", False)),
                )
                for bond in sorted(bonds, key=lambda item: item.id)
            ),
            tuple(rounded_coords),
        )
    ).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()


def _candidate_from_motif_constraints(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    """Translate rigid covalent components to satisfy semantic interactions.

    This is intentionally not an atom-level planar optimizer: each covalent
    component is treated as a block and only whole-block translations are
    applied before normal candidate ranking expands the result to atom coords.
    """
    layer_model = build_multilayer_chemical_graph(graph, atom_ids)
    constraints = tuple(
        constraint
        for constraint in layer_model.layout_constraint_graph.constraints
        if constraint.kind == LayoutConstraintKind.INTERACTION_DISTANCE
        and len(constraint.atom_ids) == 2
    )
    if not constraints:
        return None

    components = _covalent_components(atom_ids, bonds)
    atom_to_component = {
        atom_id: idx
        for idx, component in enumerate(components)
        for atom_id in component
    }
    if len(components) < 2:
        return None

    coords = dict(before)
    before_error = _interaction_constraint_error(constraints, coords)
    moved_components: set[int] = set()

    for _ in range(5):
        offsets = {idx: (0.0, 0.0, 0.0) for idx in range(len(components))}
        for constraint in constraints:
            a_id, b_id = constraint.atom_ids
            comp_a = atom_to_component.get(a_id)
            comp_b = atom_to_component.get(b_id)
            if comp_a is None or comp_b is None or comp_a == comp_b:
                continue
            if a_id not in coords or b_id not in coords:
                continue
            ax, ay = coords[a_id]
            bx, by = coords[b_id]
            dx = bx - ax
            dy = by - ay
            dist = math.hypot(dx, dy)
            if dist <= 1e-6:
                ux, uy = 1.0, 0.0
                dist = 1.0
            else:
                ux, uy = dx / dist, dy / dist
            desired = float(constraint.target_distance or target * 1.35)
            delta = dist - desired
            if abs(delta) < target * 0.04:
                continue
            step = max(-target * 0.75, min(target * 0.75, delta * 0.55))
            weight = max(0.1, float(constraint.weight or 1.0))
            move_a, move_b = _interaction_component_move_weights(
                graph,
                components[comp_a],
                components[comp_b],
                a_id,
                b_id,
            )
            if move_a:
                ox, oy, ow = offsets[comp_a]
                offsets[comp_a] = (
                    ox + ux * step * move_a * weight,
                    oy + uy * step * move_a * weight,
                    ow + weight,
                )
            if move_b:
                ox, oy, ow = offsets[comp_b]
                offsets[comp_b] = (
                    ox - ux * step * move_b * weight,
                    oy - uy * step * move_b * weight,
                    ow + weight,
                )

        changed = False
        for comp_idx, (ox, oy, weight) in offsets.items():
            if weight <= 0.0:
                continue
            offset = (ox / weight, oy / weight)
            if math.hypot(*offset) < 0.05:
                continue
            moved_components.add(comp_idx)
            changed = True
            for atom_id in components[comp_idx]:
                x, y = coords[atom_id]
                coords[atom_id] = (x + offset[0], y + offset[1])
        if not changed:
            break

    after_error = _interaction_constraint_error(constraints, coords)
    if after_error >= before_error * 0.80 or _mean_displacement(before, coords, atom_ids) < 0.25:
        return None
    return Clean2DCandidate(
        source="motif_constraints",
        coords=coords,
        message="Layout 2D ajustado por motivos e interacciones",
        metadata={
            "multilayer_clean2d": True,
            "motif_first": True,
            "interaction_constraint_count": len(constraints),
            "interaction_constraint_error": after_error,
            "interaction_constraint_error_before": before_error,
            "moved_motif_blocks": len(moved_components),
        },
    )


def _candidate_from_block_constraints(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    layer_model = build_multilayer_chemical_graph(graph, atom_ids)
    block_graph = layer_model.block_graph
    explicit_constraints = tuple(
        constraint
        for constraint in layer_model.layout_constraint_graph.constraints
        if constraint.kind == LayoutConstraintKind.INTERACTION_DISTANCE
        and len(constraint.atom_ids) == 2
    )
    block_constraints = _internal_block_layout_constraints(block_graph, before, target)
    constraints = tuple((*explicit_constraints, *block_constraints))
    quality = classify_clean2d_layout_quality(graph, atom_ids, coords=before, target_bond_length=target)
    if not block_graph.blocks or not has_intramolecular_block_layout_problem(block_graph, quality):
        return None

    components = _covalent_components(atom_ids, bonds)
    atom_component = {
        atom_id: idx
        for idx, component in enumerate(components)
        for atom_id in component
    }
    coords = dict(before)
    before_error = _block_constraint_error(block_graph, constraints, coords, bonds, target)

    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    operations: list[dict[str, object]] = []
    for _ in range(3):
        improved = False
        for edge in getattr(block_graph, "edges", ()) or ():
            if getattr(edge, "kind", None) not in {BlockEdgeKind.LINKER, BlockEdgeKind.ATTACHMENT, BlockEdgeKind.BRIDGE}:
                continue
            atom_pair = tuple(getattr(edge, "atom_ids", ()) or ())
            if len(atom_pair) != 2:
                continue
            left, right = atom_pair
            if left not in coords or right not in coords or atom_component.get(left) != atom_component.get(right):
                continue
            trial = _best_block_edge_transform(
                block_graph,
                coords,
                atom_ids,
                adjacency,
                constraints,
                bonds,
                target,
                left,
                right,
            )
            if trial is None:
                continue
            trial_coords, operation = trial
            if _block_constraint_error(block_graph, constraints, trial_coords, bonds, target) + 0.5 < _block_constraint_error(
                block_graph,
                constraints,
                coords,
                bonds,
                target,
            ):
                coords = trial_coords
                operations.append(operation)
                improved = True
        if not improved:
            break

    after_error = _block_constraint_error(block_graph, constraints, coords, bonds, target)
    if after_error >= before_error * 0.98 or _mean_displacement(before, coords, atom_ids) < 0.25:
        return None
    block_counts: dict[str, int] = {}
    for block in block_graph.blocks:
        block_counts[block.kind.value] = block_counts.get(block.kind.value, 0) + 1
    return Clean2DCandidate(
        source="block_constraints",
        coords=coords,
        message="Layout 2D ajustado por bloques intramoleculares",
        metadata={
            "multilayer_clean2d": True,
            "block_first": True,
            "block_constraint_count": len(constraints),
            "internal_block_constraint_count": len(block_constraints),
            "block_operation_count": len(operations),
            "block_kinds": block_counts,
            "block_operations": tuple(operations),
            "interaction_constraint_count": len(explicit_constraints),
            "interaction_constraint_error": _interaction_constraint_error(explicit_constraints, coords),
            "interaction_constraint_error_before": _interaction_constraint_error(explicit_constraints, before),
            "block_constraint_error": after_error,
            "block_constraint_error_before": before_error,
        },
    )


def _candidate_from_block_layout_signals(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    layer_model = build_multilayer_chemical_graph(graph, atom_ids)
    if not has_hierarchical_block_layout_signals(layer_model.block_graph):
        return None
    candidate = _candidate_from_block_constraints(graph, atom_ids, before, bonds, target)
    if candidate is None:
        block_constraints = _internal_block_layout_constraints(layer_model.block_graph, before, target)
        before_error = _block_constraint_error(layer_model.block_graph, block_constraints, before, bonds, target)
        after = _normalized_candidate_coords(before, before, bonds, target)
        candidate_metadata = {
            "multilayer_clean2d": True,
            "block_first": True,
            "block_layout_fallback": "rigid_global_normalization",
            "block_constraint_count": len(block_constraints),
            "internal_block_constraint_count": len(block_constraints),
            "block_operation_count": 0,
            "interaction_constraint_count": 0,
            "interaction_constraint_error": 0.0,
            "interaction_constraint_error_before": 0.0,
            "block_constraint_error": _block_constraint_error(layer_model.block_graph, block_constraints, after, bonds, target),
            "block_constraint_error_before": before_error,
        }
    else:
        after = _complete_coords(candidate.coords, before, atom_ids)
        candidate_metadata = candidate.metadata
    rejection = _reject_block_layout_regression(graph, atom_ids, bonds, before, after, target)
    if rejection:
        return Clean2DCandidate(
            source="block_layout",
            coords=after,
            rejected=True,
            rejection_reason=rejection,
            metadata={**candidate_metadata, "block_layout_rejected": rejection},
        )
    return Clean2DCandidate(
        source="block_layout",
        coords=after,
        message="Layout 2D jerarquico por BlockGraph",
        metadata={
            **candidate_metadata,
            "block_layout": True,
            "hierarchical_block_signals": True,
        },
    )


def _reject_block_layout_regression(
    graph: MolGraph,
    atom_ids: set[int],
    bonds: list[Bond],
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    target: float,
) -> str:
    before_crossings = _count_crossings(before, bonds)
    after_crossings = _count_crossings(after, bonds)
    if after_crossings > before_crossings:
        return "block_layout_aumenta_cruces"
    before_nonbonded = min_nonbonded_distance(before, bonds, atom_ids)
    after_nonbonded = min_nonbonded_distance(after, bonds, atom_ids)
    if before_nonbonded < math.inf and after_nonbonded + 1e-6 < min(before_nonbonded, target * 0.35):
        return "block_layout_empeora_colisiones"
    if stereo_layout_signature(graph, before, atom_ids) != stereo_layout_signature(graph, after, atom_ids):
        return "block_layout_cambia_firma_estereo"
    before_quality = classify_clean2d_layout_quality(graph, atom_ids, coords=before, target_bond_length=target)
    after_quality = classify_clean2d_layout_quality(graph, atom_ids, coords=after, target_bond_length=target)
    if after_quality.min_ring_degeneracy + 1e-6 < before_quality.min_ring_degeneracy:
        return "block_layout_colapsa_anillos"
    return ""


def _candidate_from_rdkit_isolated(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
    timeout_s: float,
) -> Clean2DCandidate | None:
    try:
        from chemuson.chemio.rdkit_safe import clean2d_isolated

        raw, error = clean2d_isolated(graph, timeout_s=timeout_s)
    except Exception as exc:
        return Clean2DCandidate(
            source="rdkit_isolated",
            coords={},
            rejected=True,
            rejection_reason=f"rdkit_aislado_no_disponible:{exc}",
        )
    if error or not raw:
        return Clean2DCandidate(
            source="rdkit_isolated",
            coords={},
            rejected=True,
            rejection_reason=f"rdkit_aislado_no_disponible:{error or 'sin_coordenadas'}",
        )
    raw = {int(atom_id): (float(x), float(y)) for atom_id, (x, y) in raw.items() if int(atom_id) in atom_ids}
    raw = _project_missing_explicit_hydrogens(graph, atom_ids, before, raw, bonds)
    coords = _normalized_candidate_coords(before, raw, bonds, target, align=True)
    return Clean2DCandidate(
        source="rdkit_isolated",
        coords=coords,
        message="Estructura 2D limpiada (RDKit aislado)",
    )


def _candidate_from_rdkit_direct(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    return Clean2DCandidate(
        source="rdkit_direct",
        coords={},
        rejected=True,
        rejection_reason="rdkit_directo_deshabilitado:usar_rdkit_aislado",
    )


def _candidate_from_internal_templates(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    coords = _internal_template_layout(graph, atom_ids, before, bonds, target)
    if not coords:
        return None
    coords = _normalized_candidate_coords(before, coords, bonds, target, align=True)
    return Clean2DCandidate(
        source="internal_templates",
        coords=coords,
        message="Estructura 2D limpiada (motor interno)",
    )


def _candidate_from_v2(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
    mode: Clean2DMode,
) -> Clean2DCandidate | None:
    try:
        params = (
            Clean2DParameters.publication(target)
            if mode == Clean2DMode.PUBLICATION
            else Clean2DParameters.quick(target)
        )
        coords = optimize_clean2d_positions(graph, atom_ids, params)
        coords = _normalized_candidate_coords(before, coords, bonds, target, align=False)
        return Clean2DCandidate(
            source="clean2d_v2",
            coords=coords,
            message="Estructura 2D limpiada (clean2d_v2)",
        )
    except Exception as exc:
        return Clean2DCandidate(
            source="clean2d_v2",
            coords={},
            rejected=True,
            rejection_reason=f"clean2d_v2_fallo:{exc}",
        )


def _candidate_from_structure_preserving_fallback(
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> Clean2DCandidate | None:
    try:
        if has_cycles(atom_ids, bonds):
            changed = structure_preserving_geometry_polish(
                before,
                bonds,
                target,
                max_displacement_per_atom=target * 8.0,
            )
        else:
            changed = structure_preserving_length_polish(
                before,
                bonds,
                target,
                max_iterations=8,
                max_displacement_per_atom=target * 8.0,
            )
            if not changed:
                changed = length_only_polish(
                    before,
                    bonds,
                    target,
                    max_iterations=36,
                    damping=0.6,
                    max_displacement_per_atom=target * 5.0,
                )
        coords = dict(before)
        coords.update(changed)
        return Clean2DCandidate(
            source="safe_fallback",
            coords=coords,
            message="Estructura 2D limpiada (fallback seguro)",
        )
    except Exception as exc:
        return Clean2DCandidate(
            source="safe_fallback",
            coords={},
            rejected=True,
            rejection_reason=f"fallback_seguro_fallo:{exc}",
        )


def _propose_candidates(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
    *,
    seed: int | None,
) -> list[Clean2DCandidate]:
    candidates: list[Clean2DCandidate] = []
    for idx, (source, coords) in enumerate(_rotatable_side_variants(atom_ids, bonds, before)):
        coords = _normalized_candidate_coords(before, coords, bonds, target, align=False)
        candidates.append(
            Clean2DCandidate(
                source=source,
                coords=coords,
                message="Estructura 2D alternativa propuesta",
                metadata={"variant_index": idx},
            )
        )

    if not _has_explicit_stereo(graph, atom_ids, bonds):
        mirrored = _mirror_coords(before, axis="x")
        mirrored = _normalized_candidate_coords(before, mirrored, bonds, target, align=False)
        candidates.append(
            Clean2DCandidate(
                source="propose_mirror",
                coords=mirrored,
                message="Estructura 2D alternativa propuesta",
            )
        )

    projection = _candidate_from_3d_projection(graph, atom_ids, before, bonds, target, seed=seed)
    if projection is not None:
        candidates.append(projection)
    return candidates


def _candidate_from_3d_projection(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
    *,
    seed: int | None,
) -> Clean2DCandidate | None:
    try:
        from chemuson.chemio.rdkit_safe import conformer_3d_isolated

        positions, _metadata, error = conformer_3d_isolated(graph, timeout_s=8.0)
        if error or not positions:
            return None
    except Exception:
        return None

    rotation_seed = int(seed if seed is not None else 0xC0FFEE)
    pitch = math.radians(18.0 + (rotation_seed % 11))
    yaw = math.radians(37.0 + (rotation_seed % 17))
    projected: dict[int, tuple[float, float]] = {}
    center = _center(before)
    cos_p = math.cos(pitch)
    sin_p = math.sin(pitch)
    cos_y = math.cos(yaw)
    sin_y = math.sin(yaw)
    for atom_id in atom_ids:
        if atom_id not in positions:
            continue
        x, y, z = positions[atom_id]
        ry_x = x * cos_y + z * sin_y
        ry_z = -x * sin_y + z * cos_y
        rp_y = y * cos_p - ry_z * sin_p
        projected[atom_id] = (ry_x, rp_y)
    if len(projected) != len(atom_ids):
        return None
    projected = _normalized_candidate_coords(before, projected, bonds, target, align=True)
    projected = _translate_to_center(projected, center)
    return Clean2DCandidate(
        source="propose_3d_projection",
        coords=projected,
        message="Estructura 2D alternativa propuesta",
    )


def _internal_template_layout(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> dict[int, tuple[float, float]]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    rings = _cycle_basis_ordered(atom_ids, bonds, max_size=8)
    rings.sort(key=lambda ring: (-len(ring), min(ring), tuple(ring)))
    out = dict(before)
    placed: set[int] = set()
    ring_centers: dict[int, tuple[float, float]] = {}
    pending_rings: list[list[int]] = []

    for ring_idx, ring in enumerate(rings):
        if not (3 <= len(ring) <= 8):
            continue
        if placed and not (set(ring) & placed):
            pending_rings.append(ring)
            continue
        coords = _place_ring_template(ring, before, out, placed, ring_centers, target, adjacency, bonds)
        if not coords:
            continue
        for atom_id, coord in coords.items():
            out[atom_id] = coord
            placed.add(atom_id)
        ring_centers[ring_idx] = _center({atom_id: out[atom_id] for atom_id in ring})

    if placed:
        deferred_rings = [set(ring) for ring in pending_rings]
        _layout_branches_from_placed_core(
            graph,
            atom_ids,
            bonds,
            adjacency,
            out,
            placed,
            target,
            deferred_rings=deferred_rings,
        )
        for ring in pending_rings:
            if not (set(ring) & placed):
                continue
            coords = _place_ring_template(ring, before, out, placed, ring_centers, target, adjacency, bonds)
            if not coords:
                continue
            for atom_id, coord in coords.items():
                out[atom_id] = coord
                placed.add(atom_id)
        _layout_branches_from_placed_core(graph, atom_ids, bonds, adjacency, out, placed, target)

    remaining = atom_ids - placed
    visited: set[int] = set()
    for start in sorted(remaining):
        if start in visited:
            continue
        component = _component_from(start, adjacency, atom_ids)
        visited.update(component)
        if component & placed:
            continue
        rebuilt = _layout_tree_component(graph, component, adjacency, before, bonds, target)
        out.update(rebuilt)
        placed.update(component)

    if not placed:
        rebuilt_all = _layout_tree_component(graph, atom_ids, adjacency, before, bonds, target)
        out.update(rebuilt_all)
    return out


def _place_ring_template(
    ring: list[int],
    before: dict[int, tuple[float, float]],
    current: dict[int, tuple[float, float]],
    placed: set[int],
    ring_centers: dict[int, tuple[float, float]],
    target: float,
    adjacency: dict[int, set[int]],
    bonds: list[Bond],
) -> dict[int, tuple[float, float]]:
    shared_edges: list[tuple[int, int, int]] = []
    n = len(ring)
    ring_set = set(ring)
    for idx, atom_id in enumerate(ring):
        nxt = ring[(idx + 1) % n]
        if atom_id in placed and nxt in placed:
            shared_edges.append((idx, atom_id, nxt))
    if shared_edges:
        idx, a1, a2 = shared_edges[0]
        ordered = ring[idx:] + ring[:idx]
        coords = _regular_ring_from_edge(ordered, current[a1], current[a2], target, sign=1)
        alt = _regular_ring_from_edge(ordered, current[a1], current[a2], target, sign=-1)
        return min(
            (coords, alt),
            key=lambda item: _ring_overlap_score(item, placed, current, shared={a1, a2}),
        )

    shared_atoms = ring_set & placed
    if 3 <= n <= 6 and len(shared_atoms) == 1:
        anchor = next(iter(shared_atoms))
        external_neighbors = sorted(
            neigh
            for neigh in adjacency.get(anchor, set())
            if neigh in placed and neigh in current and neigh not in ring_set
        )
        if external_neighbors:
            external = external_neighbors[0]
            candidates = [
                _place_ring_from_single_anchor(ring, anchor, external, current, target, sign=sign)
                for sign in (1, -1)
            ]
            candidates = [coords for coords in candidates if coords]
            if candidates:
                return min(
                    candidates,
                    key=lambda item: _single_anchor_ring_orientation_score(
                        item,
                        ring,
                        anchor,
                        external,
                        current,
                        placed,
                        bonds,
                        target,
                    ),
                )

    center = _center({atom_id: before[atom_id] for atom_id in ring if atom_id in before})
    if center == (0.0, 0.0) and not any(atom_id in before for atom_id in ring):
        center = (0.0, 0.0)
    coords = _regular_ring_at_center(ring, center, target)
    if placed:
        shared = [atom_id for atom_id in ring if atom_id in placed]
        if shared:
            anchor = shared[0]
            ax, ay = coords[anchor]
            tx, ty = current[anchor]
            coords = {atom_id: (x + tx - ax, y + ty - ay) for atom_id, (x, y) in coords.items()}
    else:
        coords = _align_to_reference({atom_id: before[atom_id] for atom_id in ring}, coords)
    return coords


def _regular_ring_at_center(
    ring: list[int],
    center: tuple[float, float],
    target: float,
) -> dict[int, tuple[float, float]]:
    n = len(ring)
    radius = target / (2.0 * math.sin(math.pi / n))
    before_angle = 0.0
    if ring:
        before_angle = math.radians(30.0 if n == 6 else -90.0)
    coords: dict[int, tuple[float, float]] = {}
    for idx, atom_id in enumerate(ring):
        angle = before_angle + 2.0 * math.pi * idx / n
        coords[atom_id] = (
            center[0] + radius * math.cos(angle),
            center[1] + radius * math.sin(angle),
        )
    return coords


def _regular_ring_from_edge(
    ordered_ring: list[int],
    p0: tuple[float, float],
    p1: tuple[float, float],
    target: float,
    *,
    sign: int,
) -> dict[int, tuple[float, float]]:
    n = len(ordered_ring)
    edge_len = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
    if edge_len <= 1e-6:
        edge_len = target
        p1 = (p0[0] + target, p0[1])
    coords: dict[int, tuple[float, float]] = {
        ordered_ring[0]: p0,
        ordered_ring[1]: p1,
    }
    angle = math.atan2(p1[1] - p0[1], p1[0] - p0[0])
    turn = float(sign) * (2.0 * math.pi / n)
    prev = p1
    for idx in range(2, n):
        angle += turn
        nxt = (prev[0] + edge_len * math.cos(angle), prev[1] + edge_len * math.sin(angle))
        coords[ordered_ring[idx]] = nxt
        prev = nxt
    return coords


def _place_ring_from_single_anchor(
    ring: list[int],
    anchor_atom_id: int,
    external_neighbor_id: int,
    current: dict[int, tuple[float, float]],
    target: float,
    *,
    sign: int,
) -> dict[int, tuple[float, float]]:
    if anchor_atom_id not in current or external_neighbor_id not in current:
        return {}
    n = len(ring)
    if n < 3:
        return {}
    anchor = current[anchor_atom_id]
    external = current[external_neighbor_id]
    vx = external[0] - anchor[0]
    vy = external[1] - anchor[1]
    norm = math.hypot(vx, vy)
    if norm <= 1e-9:
        return {}
    ux = vx / norm
    uy = vy / norm
    radius = target / (2.0 * math.sin(math.pi / n))
    center = (anchor[0] - ux * radius, anchor[1] - uy * radius)
    anchor_angle = math.atan2(anchor[1] - center[1], anchor[0] - center[0])
    if anchor_atom_id not in ring:
        return {}
    start_idx = ring.index(anchor_atom_id)
    ordered = ring[start_idx:] + ring[:start_idx]
    coords: dict[int, tuple[float, float]] = {}
    for idx, atom_id in enumerate(ordered):
        angle = anchor_angle + float(sign) * 2.0 * math.pi * idx / n
        coords[atom_id] = (
            center[0] + radius * math.cos(angle),
            center[1] + radius * math.sin(angle),
        )
    coords[anchor_atom_id] = anchor
    return coords


def _single_anchor_ring_orientation_score(
    coords: dict[int, tuple[float, float]],
    ring: list[int],
    anchor: int,
    external_neighbor: int,
    current: dict[int, tuple[float, float]],
    placed: set[int],
    bonds: list[Bond],
    target: float,
) -> tuple[float, ...]:
    if anchor not in coords or external_neighbor not in current:
        return (float("inf"),)
    ring_center = _center({atom_id: coords[atom_id] for atom_id in ring if atom_id in coords})
    anchor_pt = coords[anchor]
    external_pt = current[external_neighbor]
    center_angle = _angle_between(anchor_pt, external_pt, ring_center)
    ring_neighbors = _ring_neighbors_for_anchor(ring, anchor)
    edge_angles = [
        _angle_between(anchor_pt, external_pt, coords[neigh])
        for neigh in ring_neighbors
        if neigh in coords
    ]
    min_edge_angle = min(edge_angles) if edge_angles else 0.0
    orientation_penalty = max(0.0, 115.0 - center_angle) * 1000.0
    edge_penalty = max(0.0, 85.0 - min_edge_angle) * 1000.0
    collisions = 0
    min_dist = float("inf")
    ring_set = set(ring)
    for atom_id, coord in coords.items():
        if atom_id == anchor:
            continue
        for placed_id in placed:
            if placed_id == anchor or placed_id in ring_set or placed_id not in current:
                continue
            dist = _distance(coord, current[placed_id])
            min_dist = min(min_dist, dist)
            if dist < target * 0.55:
                collisions += 1
    crossings = _single_anchor_ring_crossings(coords, ring, current, placed, bonds)
    overlap = _ring_overlap_score(coords, placed, current, shared={anchor})
    return (
        orientation_penalty,
        edge_penalty,
        float(collisions),
        float(crossings),
        -center_angle,
        -min_edge_angle,
        -min(min_dist, target * 3.0),
        float(overlap[0]),
        float(overlap[1]),
    )


def _ring_neighbors_for_anchor(ring: list[int], anchor: int) -> list[int]:
    if anchor not in ring:
        return []
    idx = ring.index(anchor)
    return [ring[(idx - 1) % len(ring)], ring[(idx + 1) % len(ring)]]


def _single_anchor_ring_crossings(
    coords: dict[int, tuple[float, float]],
    ring: list[int],
    current: dict[int, tuple[float, float]],
    placed: set[int],
    bonds: list[Bond],
) -> int:
    crossings = 0
    ring_edges = [
        (ring[idx], ring[(idx + 1) % len(ring)])
        for idx in range(len(ring))
        if ring[idx] in coords and ring[(idx + 1) % len(ring)] in coords
    ]
    for a1, a2 in ring_edges:
        for bond in bonds:
            if bond.a1_id not in placed or bond.a2_id not in placed:
                continue
            if bond.a1_id not in current or bond.a2_id not in current:
                continue
            if {a1, a2} & {bond.a1_id, bond.a2_id}:
                continue
            if _segments_intersect(coords[a1], coords[a2], current[bond.a1_id], current[bond.a2_id]):
                crossings += 1
    return crossings


def _ring_overlap_score(
    coords: dict[int, tuple[float, float]],
    placed: set[int],
    current: dict[int, tuple[float, float]],
    *,
    shared: set[int],
) -> tuple[int, float]:
    collisions = 0
    min_dist = float("inf")
    for atom_id, coord in coords.items():
        if atom_id in shared:
            continue
        for placed_id in placed:
            if placed_id in shared:
                continue
            dist = _distance(coord, current[placed_id])
            min_dist = min(min_dist, dist)
            if dist < 12.0:
                collisions += 1
    return collisions, -min_dist


def _layout_branches_from_placed_core(
    graph: MolGraph,
    atom_ids: set[int],
    bonds: list[Bond],
    adjacency: dict[int, set[int]],
    out: dict[int, tuple[float, float]],
    placed: set[int],
    target: float,
    *,
    deferred_rings: list[set[int]] | None = None,
) -> None:
    deferred_rings = deferred_rings or []
    ring_center_by_atom = _ring_center_by_atom(placed, adjacency, out)
    queue = sorted(placed)
    parent: dict[int, int | None] = {atom_id: None for atom_id in placed}
    while queue:
        anchor = queue.pop(0)
        unplaced = sorted(
            neigh
            for neigh in adjacency.get(anchor, set())
            if neigh in atom_ids
            and neigh not in placed
            and not _defer_pending_ring_edge(anchor, neigh, deferred_rings)
        )
        if not unplaced:
            continue
        base_angles: list[float] = []
        ring_center = ring_center_by_atom.get(anchor)
        if ring_center is not None:
            ax, ay = out[anchor]
            base_angles.append(math.degrees(math.atan2(ay - ring_center[1], ax - ring_center[0])))
        parent_id = parent.get(anchor)
        if parent_id is not None and parent_id in out:
            incoming = _angle_deg(out[anchor], out[parent_id])
            geometry = _atom_geometry(graph, anchor, bonds)
            spread = 180.0 if geometry == "sp" else 120.0
            base_angles.extend([(incoming + spread) % 360.0, (incoming - spread) % 360.0])
        if not base_angles:
            base_angles.append(0.0)

        for idx, neigh in enumerate(unplaced):
            bond = _bond_between(bonds, anchor, neigh)
            length = _target_length_for_bond(bond, target)
            angle = _choose_branch_angle(anchor, neigh, base_angles, out, placed, bonds, length, target)
            out[neigh] = _endpoint(out[anchor], angle, length)
            placed.add(neigh)
            parent[neigh] = anchor
            queue.append(neigh)
            base_angles.append(angle)


def _defer_pending_ring_edge(anchor: int, neighbor: int, deferred_rings: list[set[int]]) -> bool:
    for ring in deferred_rings:
        if anchor in ring and neighbor in ring:
            return True
    return False


def _layout_tree_component(
    graph: MolGraph,
    component: set[int],
    adjacency: dict[int, set[int]],
    before: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> dict[int, tuple[float, float]]:
    if not component:
        return {}
    root = _tree_root(graph, component, adjacency, before)
    out: dict[int, tuple[float, float]] = {root: before.get(root, (0.0, 0.0))}
    placed = {root}
    parent: dict[int, int | None] = {root: None}
    queue = [root]
    root_neighbors = sorted(adjacency.get(root, set()) & component)
    root_angles = _spread_angles(len(root_neighbors), 0.0)
    while queue:
        node = queue.pop(0)
        neighbors = sorted(neigh for neigh in adjacency.get(node, set()) if neigh in component and neigh not in placed)
        if not neighbors:
            continue
        existing_angles: list[float] = []
        if parent.get(node) is not None:
            incoming = _angle_deg(out[node], out[parent[node]])  # type: ignore[index]
            geometry = _atom_geometry(graph, node, bonds)
            spread = 180.0 if geometry == "sp" else 120.0
            existing_angles.extend([(incoming + spread) % 360.0, (incoming - spread) % 360.0])
        elif node == root:
            existing_angles.extend(root_angles)
        if not existing_angles:
            existing_angles.append(0.0)
        for idx, neigh in enumerate(neighbors):
            bond = _bond_between(bonds, node, neigh)
            length = _target_length_for_bond(bond, target)
            base = existing_angles[idx % len(existing_angles)]
            angle = _choose_branch_angle(node, neigh, [base, base + 60.0, base - 60.0], out, placed, bonds, length, target)
            out[neigh] = _endpoint(out[node], angle, length)
            placed.add(neigh)
            parent[neigh] = node
            queue.append(neigh)

    before_center = _center({atom_id: before[atom_id] for atom_id in component if atom_id in before})
    return _translate_to_center(out, before_center)


def _choose_branch_angle(
    anchor: int,
    new_atom: int,
    base_angles: list[float],
    out: dict[int, tuple[float, float]],
    placed: set[int],
    bonds: list[Bond],
    length: float,
    target: float,
) -> float:
    candidates: list[float] = []
    for angle in base_angles:
        for delta in (0.0, 60.0, -60.0, 120.0, -120.0, 180.0):
            candidates.append((angle + delta) % 360.0)
    for grid in range(0, 360, 30):
        candidates.append(float(grid))

    seen: set[int] = set()
    unique: list[float] = []
    for angle in candidates:
        key = int(round(angle * 10.0))
        if key in seen:
            continue
        seen.add(key)
        unique.append(angle)

    def score(angle: float) -> tuple[int, float]:
        point = _endpoint(out[anchor], angle, length)
        collisions = 0
        min_dist = float("inf")
        for atom_id in placed:
            if atom_id == anchor:
                continue
            dist = _distance(point, out[atom_id])
            min_dist = min(min_dist, dist)
            if dist < target * 0.55:
                collisions += 1
        crossings = 0
        for bond in bonds:
            if bond.a1_id not in placed or bond.a2_id not in placed:
                continue
            if anchor in {bond.a1_id, bond.a2_id}:
                continue
            if _segments_intersect(out[anchor], point, out[bond.a1_id], out[bond.a2_id]):
                crossings += 1
        closeness = min(_angle_distance(angle, base) for base in base_angles)
        return collisions + crossings * 4, closeness - min(min_dist, target * 2.0) * 0.05

    return min(unique, key=score)


def _rotatable_side_variants(
    atom_ids: set[int],
    bonds: list[Bond],
    before: dict[int, tuple[float, float]],
) -> list[tuple[str, dict[int, tuple[float, float]]]]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    ring_core = _cycle_core_atoms(atom_ids, adjacency)
    variants: list[tuple[tuple[int, int, int], str, dict[int, tuple[float, float]]]] = []
    for bond in bonds:
        if not _is_rotatable_bond(bond, atom_ids):
            continue
        left = _component_without_edge(int(bond.a1_id), int(bond.a2_id), adjacency, atom_ids)
        right = _component_without_edge(int(bond.a2_id), int(bond.a1_id), adjacency, atom_ids)
        if not left or not right or left & right:
            continue
        side = _preferred_side_for_variant(left, right, ring_core)
        if not side:
            continue
        reflected = _reflect_side_across_bond(before, bond, side)
        if reflected:
            coords = dict(before)
            coords.update(reflected)
            variants.append(((0, len(side), int(bond.id)), "propose_reflection", coords))
        for delta in (60.0, -60.0, 120.0, -120.0):
            rotated = _rotate_side_around_anchor(before, bond, side, delta)
            if rotated:
                coords = dict(before)
                coords.update(rotated)
                variants.append(((1, len(side), int(bond.id)), "propose_rotation", coords))
    variants.sort(key=lambda item: item[0])
    return [(source, coords) for _score, source, coords in variants]


def _prepare_rdkit_mol_for_depiction(Chem, mol) -> str:
    try:
        Chem.SanitizeMol(mol)
        return "sanitize"
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL
        if hasattr(Chem.SanitizeFlags, "SANITIZE_KEKULIZE"):
            sanitize_ops = sanitize_ops ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
        return "sanitize_without_kekulize"
    except Exception:
        pass
    try:
        Chem.SetAromaticity(mol)
        return "set_aromaticity_only"
    except Exception:
        return "unsanitized"


def _safety_mode_for_candidate(
    mode: Clean2DMode,
    candidate: Clean2DCandidate,
    baseline_bad: bool,
) -> str:
    if mode == Clean2DMode.PROPOSE:
        return "conformer"
    if mode == Clean2DMode.PUBLICATION:
        return "publication"
    if baseline_bad and candidate.source != "current":
        return "conformer"
    return "quick"


def _layout_needs_rebuild(
    atom_ids: set[int],
    bonds: list[Bond],
    coords: dict[int, tuple[float, float]],
    target: float,
) -> bool:
    if not atom_ids or not bonds:
        return False
    stats = _bond_stats(coords, bonds)
    if stats["mean"] <= 1e-6:
        return True
    if stats["min"] < target * 0.45 or stats["max"] > target * 1.85:
        return True
    if min_nonbonded_distance(coords, bonds, atom_ids) < target * 0.35:
        return True
    if _count_crossings(coords, bonds) > 0:
        return True
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    small_ring_centers = _small_ring_atoms(atom_ids, bonds)
    for center, neighbors_set in adjacency.items():
        if center in small_ring_centers:
            continue
        neighbors = [neigh for neigh in neighbors_set if neigh in coords]
        if len(neighbors) < 2 or center not in coords:
            continue
        if len(neighbors) >= 4:
            continue
        center_bonds = [bond for bond in bonds if center in {bond.a1_id, bond.a2_id}]
        has_multiple = any(
            int(getattr(bond, "order", 1) or 1) >= 2 or bool(getattr(bond, "is_aromatic", False))
            for bond in center_bonds
        )
        for idx, left in enumerate(neighbors):
            for right in neighbors[idx + 1 :]:
                angle = _angle_between(coords[center], coords[left], coords[right])
                if not has_multiple and angle > 155.0:
                    return True
                if has_multiple and (angle < 75.0 or angle > 170.0):
                    return True
    if has_cycles(atom_ids, bonds):
        rings = _cycle_basis_ordered(atom_ids, bonds, max_size=8)
        for ring in rings:
            if ring_degeneracy_score(coords, set(ring)) < 0.12:
                return True
    return False


def _visual_quality_score(
    graph: MolGraph,
    atom_ids: set[int],
    bonds: list[Bond],
    before: dict[int, tuple[float, float]],
    coords: dict[int, tuple[float, float]],
    target: float,
    mode: Clean2DMode,
) -> float:
    score = 0.0
    lengths = []
    for bond in bonds:
        if bond.a1_id not in coords or bond.a2_id not in coords:
            score += target * 1000.0
            continue
        desired = _target_length_for_bond(bond, target)
        length = _distance(coords[bond.a1_id], coords[bond.a2_id])
        lengths.append(length)
        score += ((length - desired) / target) ** 2 * 120.0
    score += _angle_quality_penalty(graph, atom_ids, bonds, coords, target)
    crossings = _count_crossings(coords, bonds)
    score += crossings * target * 60.0
    nonbonded = min_nonbonded_distance(coords, bonds, atom_ids)
    if nonbonded < target * 0.55:
        score += (target * 0.55 - nonbonded) * 25.0
    for ring in _cycle_basis_ordered(atom_ids, bonds, max_size=8):
        degeneracy = ring_degeneracy_score(coords, set(ring))
        if degeneracy < 0.35:
            score += (0.35 - degeneracy) * target * 120.0
    exocyclic_stats = _exocyclic_ring_orientation_stats(atom_ids, bonds, coords)
    score += exocyclic_stats["bad_count"] * target * 500.0
    if exocyclic_stats["min_center_angle"] < math.inf:
        score += max(0.0, 125.0 - exocyclic_stats["min_center_angle"]) * target * 2.0
    if exocyclic_stats["min_edge_angle"] < math.inf:
        score += max(0.0, 90.0 - exocyclic_stats["min_edge_angle"]) * target * 2.0
    bbox_ratio = _bbox_ratio(before, coords, atom_ids)
    if bbox_ratio > 2.5:
        score += (bbox_ratio - 2.5) * target * 5.0
    elif bbox_ratio < 0.4:
        score += (0.4 - bbox_ratio) * target * 5.0
    return score


def _angle_quality_penalty(
    graph: MolGraph,
    atom_ids: set[int],
    bonds: list[Bond],
    coords: dict[int, tuple[float, float]],
    target: float,
) -> float:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    penalty = 0.0
    for center in sorted(atom_ids):
        neighbors = [neigh for neigh in adjacency.get(center, set()) if neigh in coords]
        if len(neighbors) < 2 or center not in coords:
            continue
        geometry = _atom_geometry(graph, center, bonds)
        preferred = 180.0 if geometry == "sp" else 120.0
        if geometry == "sp3":
            preferred = 109.5
        if len(neighbors) == 2 and any(_bond_between(bonds, center, neigh) and bool(getattr(_bond_between(bonds, center, neigh), "is_aromatic", False)) for neigh in neighbors):
            preferred = 120.0
        for idx, left in enumerate(neighbors):
            for right in neighbors[idx + 1 :]:
                angle = _angle_between(coords[center], coords[left], coords[right])
                diff = min(abs(angle - preferred), abs(angle - (360.0 - preferred)))
                if geometry == "sp3":
                    diff = min(diff, abs(angle - 120.0))
                penalty += (diff / 30.0) ** 2 * 8.0
    return penalty


def _strict_angle_deviation_stats(
    graph: MolGraph,
    atom_ids: set[int],
    bonds: list[Bond],
    coords: dict[int, tuple[float, float]],
) -> dict[str, float]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    small_ring_triples = _small_ring_angle_triples(atom_ids, bonds)
    small_ring_centers = _small_ring_atoms(atom_ids, bonds)
    deviations: list[float] = []
    penalty = 0.0
    severe_count = 0
    for center in sorted(atom_ids):
        if center in small_ring_centers:
            continue
        neighbors = [neigh for neigh in adjacency.get(center, set()) if neigh in coords]
        if len(neighbors) < 2 or center not in coords:
            continue
        if len(neighbors) >= 4:
            continue
        geometry = _atom_geometry(graph, center, bonds)
        preferred = _preferred_angle_for_center(graph, center, neighbors, bonds, geometry)
        for idx, left in enumerate(neighbors):
            for right in neighbors[idx + 1 :]:
                if frozenset((left, center, right)) in small_ring_triples:
                    continue
                angle = _angle_between(coords[center], coords[left], coords[right])
                diff = min(abs(angle - preferred), abs(angle - (360.0 - preferred)))
                if geometry == "sp3":
                    diff = min(diff, abs(angle - 120.0))
                deviations.append(diff)
                penalty += (diff / 30.0) ** 2 * 8.0
                if diff > 34.0:
                    severe_count += 1
    if not deviations:
        return {"rms": 0.0, "max": 0.0, "penalty": 0.0, "severe_count": 0.0}
    return {
        "rms": math.sqrt(sum(diff * diff for diff in deviations) / len(deviations)),
        "max": max(deviations),
        "penalty": penalty,
        "severe_count": float(severe_count),
    }


def _exocyclic_ring_orientation_stats(
    atom_ids: set[int],
    bonds: list[Bond],
    coords: dict[int, tuple[float, float]],
) -> dict[str, float]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    min_center_angle = math.inf
    min_edge_angle = math.inf
    bad_count = 0
    for ring in _cycle_basis_ordered(atom_ids, bonds, max_size=6):
        ring_set = set(ring)
        anchors = [
            atom_id
            for atom_id in ring
            if any(neigh not in ring_set for neigh in adjacency.get(atom_id, set()))
        ]
        if len(anchors) != 1:
            continue
        anchor = anchors[0]
        if anchor not in coords:
            continue
        external_neighbors = [
            neigh
            for neigh in sorted(adjacency.get(anchor, set()) - ring_set)
            if neigh in coords
        ]
        if not external_neighbors:
            continue
        external = external_neighbors[0]
        ring_coords = {atom_id: coords[atom_id] for atom_id in ring if atom_id in coords}
        if len(ring_coords) != len(ring):
            continue
        ring_center = _center(ring_coords)
        center_angle = _angle_between(coords[anchor], coords[external], ring_center)
        edge_angles = [
            _angle_between(coords[anchor], coords[external], coords[neigh])
            for neigh in _ring_neighbors_for_anchor(ring, anchor)
            if neigh in coords
        ]
        edge_angle = min(edge_angles) if edge_angles else math.inf
        min_center_angle = min(min_center_angle, center_angle)
        min_edge_angle = min(min_edge_angle, edge_angle)
        if center_angle < 110.0 or edge_angle < 85.0:
            bad_count += 1
    return {
        "min_center_angle": min_center_angle,
        "min_edge_angle": min_edge_angle,
        "bad_count": float(bad_count),
    }


def _preferred_angle_for_center(
    graph: MolGraph,
    center: int,
    neighbors: list[int],
    bonds: list[Bond],
    geometry: str,
) -> float:
    preferred = 180.0 if geometry == "sp" else 120.0
    if geometry == "sp3":
        preferred = 109.5
    if len(neighbors) == 2 and any(
        _bond_between(bonds, center, neigh)
        and bool(getattr(_bond_between(bonds, center, neigh), "is_aromatic", False))
        for neigh in neighbors
    ):
        preferred = 120.0
    return preferred


def _small_ring_angle_triples(atom_ids: set[int], bonds: list[Bond]) -> set[frozenset[int]]:
    triples: set[frozenset[int]] = set()
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    for ring in _cycle_basis_ordered(atom_ids, bonds, max_size=4):
        ring_set = set(ring)
        for center in ring:
            ring_neighbors = sorted(adjacency.get(center, set()) & ring_set)
            for idx, left in enumerate(ring_neighbors):
                for right in ring_neighbors[idx + 1 :]:
                    triples.add(frozenset((left, center, right)))
    return triples


def _small_ring_atoms(atom_ids: set[int], bonds: list[Bond]) -> set[int]:
    atoms: set[int] = set()
    for ring in _cycle_basis_ordered(atom_ids, bonds, max_size=4):
        atoms.update(ring)
    return atoms


def _deduplicate_candidates(
    graph: MolGraph,
    atom_ids: set[int],
    candidates: list[Clean2DCandidate],
    *,
    avoid_hashes: set[str] | None,
) -> list[Clean2DCandidate]:
    seen: set[str] = set()
    deduped: list[Clean2DCandidate] = []
    avoid_hashes = set(avoid_hashes or set())
    for candidate in candidates:
        if candidate.rejected and not candidate.coords:
            deduped.append(candidate)
            continue
        geometry_hash = clean2d_geometry_hash(graph, candidate.coords, atom_ids)
        if geometry_hash in seen and candidate.source != "current":
            continue
        seen.add(geometry_hash)
        deduped.append(_replace_candidate(candidate, geometry_hash=geometry_hash))
    return deduped


def _replace_candidate(candidate: Clean2DCandidate, **updates: Any) -> Clean2DCandidate:
    data = {
        "source": candidate.source,
        "coords": candidate.coords,
        "message": candidate.message,
        "score": candidate.score,
        "novelty": candidate.novelty,
        "report": candidate.report,
        "rejected": candidate.rejected,
        "rejection_reason": candidate.rejection_reason,
        "geometry_hash": candidate.geometry_hash,
        "metadata": candidate.metadata,
    }
    data.update(updates)
    return Clean2DCandidate(**data)


def _source_priority(source: str) -> int:
    priorities = {
        "current": 0,
        "block_layout": 1,
        "block_constraints": 2,
        "rdkit_isolated": 3,
        "rdkit_direct": 4,
        "internal_templates": 5,
        "clean2d_v2": 6,
        "motif_constraints": 7,
        "safe_fallback": 8,
        "local_graph": 9,
        "propose_reflection": 1,
        "propose_rotation": 2,
        "propose_3d_projection": 3,
        "propose_mirror": 4,
    }
    return priorities.get(source, 99)


def _propose_source_priority(source: str) -> int:
    if source.startswith("propose_"):
        return 0
    return _source_priority(source) + 10


def _normalize_atom_ids(graph: MolGraph, atom_ids: Iterable[int] | None) -> set[int]:
    if atom_ids is None:
        return set(graph.atoms.keys())
    return {int(atom_id) for atom_id in atom_ids if int(atom_id) in graph.atoms}


def _covalent_components(atom_ids: set[int], bonds: list[Bond]) -> list[frozenset[int]]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    seen: set[int] = set()
    components: list[frozenset[int]] = []
    for atom_id in sorted(atom_ids):
        if atom_id in seen:
            continue
        component = _component_from(atom_id, adjacency, atom_ids)
        seen.update(component)
        components.append(frozenset(component))
    return components


def _interaction_component_move_weights(
    graph: MolGraph,
    comp_a: frozenset[int],
    comp_b: frozenset[int],
    a_id: int,
    b_id: int,
) -> tuple[float, float]:
    a_center = bool(getattr(graph.atoms.get(a_id), "is_coordination_center", False))
    b_center = bool(getattr(graph.atoms.get(b_id), "is_coordination_center", False))
    if a_center and not b_center:
        return (0.0, 1.0)
    if b_center and not a_center:
        return (1.0, 0.0)
    if len(comp_a) > len(comp_b):
        return (0.0, 1.0)
    if len(comp_b) > len(comp_a):
        return (1.0, 0.0)
    return (0.5, 0.5)


def _interaction_constraint_count(constraints: Iterable[object]) -> int:
    return sum(
        1
        for constraint in constraints
        if getattr(constraint, "kind", None) == LayoutConstraintKind.INTERACTION_DISTANCE
        and len(getattr(constraint, "atom_ids", ())) == 2
    )


def _interaction_constraint_error(
    constraints: Iterable[object],
    coords: dict[int, tuple[float, float]],
) -> float:
    total = 0.0
    weight_total = 0.0
    for constraint in constraints:
        if getattr(constraint, "kind", None) != LayoutConstraintKind.INTERACTION_DISTANCE:
            continue
        atom_pair = tuple(getattr(constraint, "atom_ids", ()))
        if len(atom_pair) != 2 or atom_pair[0] not in coords or atom_pair[1] not in coords:
            continue
        target = getattr(constraint, "target_distance", None)
        if target is None:
            continue
        ax, ay = coords[atom_pair[0]]
        bx, by = coords[atom_pair[1]]
        weight = max(0.1, float(getattr(constraint, "weight", 1.0) or 1.0))
        total += abs(math.hypot(bx - ax, by - ay) - float(target)) * weight
        weight_total += weight
    if weight_total <= 0.0:
        return 0.0
    return total / weight_total


def _internal_block_layout_constraints(
    block_graph: object,
    coords: dict[int, tuple[float, float]],
    target: float,
) -> tuple[LayoutConstraint, ...]:
    constraints: list[LayoutConstraint] = []
    blocks = tuple(getattr(block_graph, "blocks", ()) or ())
    aromatic = [block for block in blocks if getattr(block, "kind", None) == BlockKind.AROMATIC_RING]
    macro_like = [
        block
        for block in blocks
        if getattr(block, "kind", None) in {BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
    ]
    cavities = [block for block in blocks if getattr(block, "kind", None) == BlockKind.INTERNAL_CAVITY]
    substituents = [block for block in blocks if getattr(block, "kind", None) == BlockKind.TERMINAL_SUBSTITUENT]
    stereo = [block for block in blocks if getattr(block, "kind", None) == BlockKind.STEREO_CENTER]

    for block in blocks:
        if getattr(block, "kind", None) in {
            BlockKind.AROMATIC_RING,
            BlockKind.FUSED_SYSTEM,
            BlockKind.MACROCYCLE,
            BlockKind.CYCLOPHANE,
        }:
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.MOTIF_RIGID,
                    atom_ids=tuple(sorted(getattr(block, "atom_ids", ()) or ())),
                    weight=4.0,
                    source=f"block_rigid:{getattr(block, 'id', 0)}",
                )
            )

    for idx, left in enumerate(aromatic):
        for right in aromatic[idx + 1 :]:
            atom_pair = _nearest_atoms_between_blocks(left, right, coords)
            if atom_pair is None:
                continue
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.INTERACTION_DISTANCE,
                    atom_ids=atom_pair,
                    target_distance=target * 1.55,
                    weight=0.7,
                    source="block_ring_centroid_separation",
                    metadata={"internal_block_constraint": "ring_centroid_separation"},
                )
            )

    for block in macro_like:
        atoms = tuple(sorted(getattr(block, "atom_ids", ()) or ()))
        if len(atoms) < 6:
            continue
        constraints.append(
            LayoutConstraint(
                id=len(constraints) + 1,
                kind=LayoutConstraintKind.MOTIF_CONTAINMENT,
                atom_ids=atoms,
                target_distance=target * min(4.0, max(2.0, len(atoms) / math.pi / 2.0)),
                weight=2.2,
                source=f"block_open_macrocycle:{getattr(block, 'id', 0)}",
                metadata={"internal_block_constraint": "open_macrocycle"},
            )
        )

    for block in cavities:
        atoms = tuple(sorted(getattr(block, "atom_ids", ()) or ()))
        if len(atoms) < 4:
            continue
        constraints.append(
            LayoutConstraint(
                id=len(constraints) + 1,
                kind=LayoutConstraintKind.MOTIF_CONTAINMENT,
                atom_ids=atoms,
                target_distance=target * 1.8,
                weight=2.0,
                source=f"block_cavity:{getattr(block, 'id', 0)}",
                metadata={"internal_block_constraint": "avoid_cavity_collapse"},
            )
        )

    for block in substituents:
        atoms = tuple(sorted(getattr(block, "atom_ids", ()) or ()))
        anchors = tuple(getattr(block, "anchor_atom_ids", ()) or ())
        external = tuple(atom_id for atom_id in atoms if atom_id not in anchors)
        if anchors and external:
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.STEREO_ORIENTATION,
                    atom_ids=(anchors[0], external[-1]),
                    target_angle_deg=120.0,
                    weight=1.0,
                    source=f"block_substituent_outward:{getattr(block, 'id', 0)}",
                    metadata={"internal_block_constraint": "substituent_outward"},
                )
            )

    for block in stereo:
        atoms = tuple(sorted(getattr(block, "atom_ids", ()) or ()))
        if atoms:
            constraints.append(
                LayoutConstraint(
                    id=len(constraints) + 1,
                    kind=LayoutConstraintKind.STEREO_ORIENTATION,
                    atom_ids=(atoms[0],),
                    weight=3.0,
                    source=f"block_stereo:{getattr(block, 'id', 0)}",
                    metadata={"internal_block_constraint": "preserve_stereo_orientation"},
                )
            )

    return tuple(constraints)


def _block_constraint_error(
    block_graph: object,
    constraints: Iterable[object],
    coords: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> float:
    total = _interaction_constraint_error(constraints, coords)
    for constraint in constraints:
        kind = getattr(constraint, "kind", None)
        atom_ids = tuple(getattr(constraint, "atom_ids", ()) or ())
        weight = max(0.1, float(getattr(constraint, "weight", 1.0) or 1.0))
        if kind == LayoutConstraintKind.MOTIF_RIGID:
            total += _rigid_block_error(atom_ids, coords, target) * weight
        elif kind == LayoutConstraintKind.MOTIF_CONTAINMENT:
            desired = float(getattr(constraint, "target_distance", None) or target * 1.5)
            total += max(0.0, desired - _block_radius(atom_ids, coords)) * weight
        elif kind == LayoutConstraintKind.STEREO_ORIENTATION:
            meta = getattr(constraint, "metadata", {}) or {}
            label = str(meta.get("internal_block_constraint", ""))
            if label == "substituent_outward" and len(atom_ids) == 2:
                total += _substituent_inward_error(block_graph, atom_ids, coords, target) * weight
    total += _linker_crossing_error(block_graph, coords, bonds, target)
    return total


def _nearest_atoms_between_blocks(left: object, right: object, coords: dict[int, tuple[float, float]]) -> tuple[int, int] | None:
    best: tuple[float, int, int] | None = None
    for a_id in getattr(left, "atom_ids", ()) or ():
        if a_id not in coords:
            continue
        for b_id in getattr(right, "atom_ids", ()) or ():
            if b_id not in coords:
                continue
            dist = _distance(coords[a_id], coords[b_id])
            if best is None or dist < best[0]:
                best = (dist, int(a_id), int(b_id))
    if best is None:
        return None
    return (best[1], best[2])


def _rigid_block_error(atom_ids: tuple[int, ...], coords: dict[int, tuple[float, float]], target: float) -> float:
    if len(atom_ids) < 3:
        return 0.0
    radius = _block_radius(atom_ids, coords)
    return max(0.0, target * 0.45 - radius)


def _block_radius(atom_ids: tuple[int, ...], coords: dict[int, tuple[float, float]]) -> float:
    present = {atom_id: coords[atom_id] for atom_id in atom_ids if atom_id in coords}
    if not present:
        return 0.0
    center = _center(present)
    return min((_distance(center, coord) for coord in present.values()), default=0.0)


def _substituent_inward_error(
    block_graph: object,
    atom_ids: tuple[int, ...],
    coords: dict[int, tuple[float, float]],
    target: float,
) -> float:
    if len(atom_ids) != 2 or atom_ids[0] not in coords or atom_ids[1] not in coords:
        return 0.0
    anchor, external = atom_ids
    parent_blocks = [
        block
        for block in getattr(block_graph, "blocks", ()) or ()
        if anchor in getattr(block, "atom_ids", frozenset())
        and getattr(block, "kind", None) in {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE}
    ]
    if not parent_blocks:
        return 0.0
    parent_center = _center({atom_id: coords[atom_id] for atom_id in parent_blocks[0].atom_ids if atom_id in coords})
    before = _distance(parent_center, coords[anchor])
    after = _distance(parent_center, coords[external])
    return max(0.0, before + target * 0.2 - after)


def _linker_crossing_error(
    block_graph: object,
    coords: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> float:
    linker_atoms = set()
    for block in getattr(block_graph, "blocks", ()) or ():
        if getattr(block, "kind", None) in {BlockKind.LINKER, BlockKind.INTRAMOLECULAR_BRIDGE}:
            linker_atoms.update(getattr(block, "atom_ids", ()) or ())
    if not linker_atoms:
        return 0.0
    linker_bonds = [bond for bond in bonds if bond.a1_id in linker_atoms or bond.a2_id in linker_atoms]
    crossings = _count_crossings(coords, linker_bonds)
    return crossings * target * 0.6


def _block_ids_for_atom(block_graph: object, atom_id: int) -> set[int]:
    mapping = getattr(block_graph, "atom_to_block_ids", {}) or {}
    return {int(block_id) for block_id in mapping.get(atom_id, ())}


def _best_internal_block_transform(
    coords: dict[int, tuple[float, float]],
    atom_ids: set[int],
    adjacency: dict[int, set[int]],
    constraints: tuple[object, ...],
    a_id: int,
    b_id: int,
    target_distance: float,
) -> tuple[dict[int, tuple[float, float]], dict[str, object]] | None:
    path = _shortest_path(adjacency, a_id, b_id, atom_ids)
    if path is None or len(path) < 4:
        return None
    current_error = _interaction_constraint_error(constraints, coords)
    current_distance_error = abs(_distance(coords[a_id], coords[b_id]) - target_distance)
    best_coords: dict[int, tuple[float, float]] | None = None
    best_operation: dict[str, object] | None = None
    best_error = current_error + current_distance_error * 0.2
    for left, right in zip(path[:-1], path[1:]):
        side = _component_without_edge(right, left, adjacency, atom_ids)
        if b_id not in side or a_id in side or len(side) >= len(atom_ids) - 1:
            continue
        pivot = coords.get(right)
        anchor = coords.get(left)
        if pivot is None or anchor is None:
            continue
        for angle in (-120.0, -90.0, -60.0, 60.0, 90.0, 120.0, 180.0):
            trial = _rotate_atom_subset(coords, side, pivot, math.radians(angle))
            error = _interaction_constraint_error(constraints, trial)
            distance_error = abs(_distance(trial[a_id], trial[b_id]) - target_distance)
            ranked_error = error + distance_error * 0.2
            if ranked_error + 0.5 < best_error:
                best_error = ranked_error
                best_coords = trial
                best_operation = {
                    "operation": "rotate",
                    "angle_deg": angle,
                    "pivot_atom_id": right,
                    "anchor_atom_id": left,
                    "moved_atom_count": len(side),
                }
        trial = _reflect_atom_subset(coords, side, anchor, pivot)
        error = _interaction_constraint_error(constraints, trial)
        distance_error = abs(_distance(trial[a_id], trial[b_id]) - target_distance)
        ranked_error = error + distance_error * 0.2
        if ranked_error + 0.5 < best_error:
            best_error = ranked_error
            best_coords = trial
            best_operation = {
                "operation": "reflect",
                "pivot_atom_id": right,
                "anchor_atom_id": left,
                "moved_atom_count": len(side),
            }
    if best_coords is None or best_operation is None:
        return None
    return best_coords, best_operation


def _best_block_edge_transform(
    block_graph: object,
    coords: dict[int, tuple[float, float]],
    atom_ids: set[int],
    adjacency: dict[int, set[int]],
    constraints: tuple[object, ...],
    bonds: list[Bond],
    target: float,
    left: int,
    right: int,
) -> tuple[dict[int, tuple[float, float]], dict[str, object]] | None:
    if left not in coords or right not in coords:
        return None
    side = _component_without_edge(right, left, adjacency, atom_ids)
    if not side or left in side or len(side) >= len(atom_ids):
        return None
    current_error = _block_constraint_error(block_graph, constraints, coords, bonds, target)
    best_coords: dict[int, tuple[float, float]] | None = None
    best_operation: dict[str, object] | None = None
    best_error = current_error
    anchor = coords[left]
    pivot = coords[right]

    for angle in (-150.0, -120.0, -90.0, -60.0, 60.0, 90.0, 120.0, 150.0, 180.0):
        trial = _rotate_atom_subset(coords, side, pivot, math.radians(angle))
        error = _block_constraint_error(block_graph, constraints, trial, bonds, target)
        if error + 0.5 < best_error:
            best_error = error
            best_coords = trial
            best_operation = {
                "operation": "rotate_block_edge",
                "angle_deg": angle,
                "anchor_atom_id": left,
                "pivot_atom_id": right,
                "moved_atom_count": len(side),
            }

    trial = _reflect_atom_subset(coords, side, anchor, pivot)
    error = _block_constraint_error(block_graph, constraints, trial, bonds, target)
    if error + 0.5 < best_error:
        best_coords = trial
        best_operation = {
            "operation": "reflect_block_edge",
            "anchor_atom_id": left,
            "pivot_atom_id": right,
            "moved_atom_count": len(side),
        }

    if best_coords is None or best_operation is None:
        return None
    return best_coords, best_operation


def _shortest_path(
    adjacency: dict[int, set[int]],
    start: int,
    end: int,
    allowed: set[int],
) -> list[int] | None:
    queue: list[tuple[int, list[int]]] = [(start, [start])]
    seen: set[int] = set()
    while queue:
        atom_id, path = queue.pop(0)
        if atom_id == end:
            return path
        if atom_id in seen:
            continue
        seen.add(atom_id)
        for neighbor in sorted(adjacency.get(atom_id, set())):
            if neighbor in allowed and neighbor not in path:
                queue.append((neighbor, path + [neighbor]))
    return None


def _rotate_atom_subset(
    coords: dict[int, tuple[float, float]],
    atom_ids: set[int],
    pivot: tuple[float, float],
    angle_rad: float,
) -> dict[int, tuple[float, float]]:
    out = dict(coords)
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    px, py = pivot
    for atom_id in atom_ids:
        if atom_id not in out:
            continue
        x, y = out[atom_id]
        dx = x - px
        dy = y - py
        out[atom_id] = (px + dx * cos_a - dy * sin_a, py + dx * sin_a + dy * cos_a)
    return out


def _reflect_atom_subset(
    coords: dict[int, tuple[float, float]],
    atom_ids: set[int],
    anchor: tuple[float, float],
    pivot: tuple[float, float],
) -> dict[int, tuple[float, float]]:
    out = dict(coords)
    ax, ay = anchor
    px, py = pivot
    ux = px - ax
    uy = py - ay
    norm = math.hypot(ux, uy)
    if norm <= 1e-6:
        return out
    ux /= norm
    uy /= norm
    for atom_id in atom_ids:
        if atom_id not in out:
            continue
        x, y = out[atom_id]
        dx = x - ax
        dy = y - ay
        parallel = dx * ux + dy * uy
        proj_x = ax + parallel * ux
        proj_y = ay + parallel * uy
        out[atom_id] = (2.0 * proj_x - x, 2.0 * proj_y - y)
    return out


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
        if bond.a1_id in atom_ids
        and bond.a2_id in atom_ids
        and bond_affects_valence(bond)
    ]


def _coerce_mode(mode: Clean2DMode | str) -> Clean2DMode:
    if isinstance(mode, Clean2DMode):
        return mode
    text = str(mode or "quick").strip().lower()
    if text in {"publication", "publish"}:
        return Clean2DMode.PUBLICATION
    if text in {"propose", "conformer", "alternative"}:
        return Clean2DMode.PROPOSE
    return Clean2DMode.QUICK


def _complete_coords(
    coords: dict[int, tuple[float, float]],
    before: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> dict[int, tuple[float, float]]:
    return {
        atom_id: (
            (float(coords[atom_id][0]), float(coords[atom_id][1]))
            if atom_id in coords
            else before[atom_id]
        )
        for atom_id in atom_ids
        if atom_id in before
    }


def _normalized_candidate_coords(
    before: dict[int, tuple[float, float]],
    coords: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
    *,
    align: bool = False,
) -> dict[int, tuple[float, float]]:
    if not coords:
        return {}
    out = {
        int(atom_id): (float(coord[0]), float(coord[1]))
        for atom_id, coord in coords.items()
        if coord is not None and len(coord) >= 2
    }
    if align:
        out = _align_to_reference(before, out)
    out = _rescale_to_bond_length(out, bonds, target)
    return _translate_to_center(out, _center(before))


def _project_missing_explicit_hydrogens(
    graph: MolGraph,
    atom_ids: set[int],
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    bonds: list[Bond],
) -> dict[int, tuple[float, float]]:
    projected = dict(after)
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    for atom_id in sorted(atom_ids):
        atom = graph.atoms.get(atom_id)
        if atom is None or atom.element != "H" or atom_id in projected:
            continue
        if atom_id not in before:
            continue
        candidates = []
        for anchor in adjacency.get(atom_id, set()):
            if anchor not in before or anchor not in projected:
                continue
            dx = before[atom_id][0] - before[anchor][0]
            dy = before[atom_id][1] - before[anchor][1]
            candidates.append((projected[anchor][0] + dx, projected[anchor][1] + dy))
        if candidates:
            projected[atom_id] = _center({idx: coord for idx, coord in enumerate(candidates)})
    return projected


def _adjacency_for_bonds(atom_ids: set[int], bonds: list[Bond]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for bond in bonds:
        if bond.a1_id in atom_ids and bond.a2_id in atom_ids:
            adjacency.setdefault(bond.a1_id, set()).add(bond.a2_id)
            adjacency.setdefault(bond.a2_id, set()).add(bond.a1_id)
    return adjacency


def _cycle_basis_ordered(
    atom_ids: set[int],
    bonds: list[Bond],
    *,
    max_size: int,
) -> list[list[int]]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    rings: list[list[int]] = []
    seen: set[tuple[int, ...]] = set()
    for bond in sorted(bonds, key=lambda item: item.id):
        path = _shortest_path_without_edge(adjacency, bond.a1_id, bond.a2_id)
        if path is None:
            continue
        if not (3 <= len(path) <= max_size):
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
) -> list[int] | None:
    blocked = (min(start, end), max(start, end))
    queue: list[tuple[int, list[int]]] = [(start, [start])]
    visited: set[int] = set()
    while queue:
        node, path = queue.pop(0)
        if node == end and len(path) >= 3:
            return path
        if node in visited:
            continue
        visited.add(node)
        for neigh in sorted(adjacency.get(node, set())):
            edge = (min(node, neigh), max(node, neigh))
            if edge == blocked:
                continue
            if neigh in path:
                continue
            queue.append((neigh, path + [neigh]))
    return None


def _components_for_snapshot(atom_ids: set[int], bonds: list[Bond]) -> tuple[tuple[int, ...], ...]:
    adjacency = _adjacency_for_bonds(atom_ids, bonds)
    components: list[tuple[int, ...]] = []
    visited: set[int] = set()
    for start in sorted(atom_ids):
        if start in visited:
            continue
        comp = _component_from(start, adjacency, atom_ids)
        visited.update(comp)
        components.append(tuple(sorted(comp)))
    return tuple(sorted(components))


def _component_from(start: int, adjacency: dict[int, set[int]], allowed: set[int]) -> set[int]:
    visited: set[int] = set()
    stack = [start]
    while stack:
        node = stack.pop()
        if node in visited or node not in allowed:
            continue
        visited.add(node)
        stack.extend(sorted(adjacency.get(node, set()) - visited))
    return visited


def _component_without_edge(
    start: int,
    blocked_neighbor: int,
    adjacency: dict[int, set[int]],
    atom_ids: set[int],
) -> set[int]:
    visited: set[int] = set()
    stack = [start]
    blocked = {start, blocked_neighbor}
    while stack:
        node = stack.pop()
        if node in visited or node not in atom_ids:
            continue
        visited.add(node)
        for neigh in adjacency.get(node, set()):
            if {node, neigh} == blocked:
                continue
            if neigh not in visited:
                stack.append(neigh)
    return visited


def _cycle_core_atoms(atom_ids: set[int], adjacency: dict[int, set[int]]) -> set[int]:
    core = set(atom_ids)
    degrees = {atom_id: len(adjacency.get(atom_id, set()) & core) for atom_id in core}
    queue = [atom_id for atom_id, degree in degrees.items() if degree <= 1]
    while queue:
        atom_id = queue.pop()
        if atom_id not in core:
            continue
        core.remove(atom_id)
        for neigh in adjacency.get(atom_id, set()):
            if neigh not in core:
                continue
            degrees[neigh] = degrees.get(neigh, 0) - 1
            if degrees[neigh] <= 1:
                queue.append(neigh)
    return core


def _is_rotatable_bond(bond: Bond, atom_ids: set[int]) -> bool:
    if bond.a1_id not in atom_ids or bond.a2_id not in atom_ids:
        return False
    if bool(getattr(bond, "is_aromatic", False)):
        return False
    if int(getattr(bond, "order", 1) or 1) != 1:
        return False
    return bond_affects_valence(bond)


def _preferred_side_for_variant(left: set[int], right: set[int], ring_core: set[int]) -> set[int]:
    left_ring = len(left & ring_core)
    right_ring = len(right & ring_core)
    if left_ring != right_ring:
        return left if left_ring < right_ring else right
    if len(left) != len(right):
        return left if len(left) < len(right) else right
    return left if min(left) <= min(right) else right


def _reflect_side_across_bond(
    before: dict[int, tuple[float, float]],
    bond: Bond,
    side: set[int],
) -> dict[int, tuple[float, float]]:
    a1 = int(bond.a1_id)
    a2 = int(bond.a2_id)
    if a1 not in before or a2 not in before:
        return {}
    x1, y1 = before[a1]
    x2, y2 = before[a2]
    dx = x2 - x1
    dy = y2 - y1
    denom = dx * dx + dy * dy
    if denom <= 1e-9:
        return {}
    reflected = {}
    for atom_id in side:
        if atom_id not in before:
            continue
        px, py = before[atom_id]
        rel_x = px - x1
        rel_y = py - y1
        projection = (rel_x * dx + rel_y * dy) / denom
        foot_x = x1 + projection * dx
        foot_y = y1 + projection * dy
        reflected[atom_id] = (2.0 * foot_x - px, 2.0 * foot_y - py)
    return reflected


def _rotate_side_around_anchor(
    before: dict[int, tuple[float, float]],
    bond: Bond,
    side: set[int],
    delta_deg: float,
) -> dict[int, tuple[float, float]]:
    a1 = int(bond.a1_id)
    a2 = int(bond.a2_id)
    if a1 not in before or a2 not in before:
        return {}
    anchor = a1 if a1 not in side else a2
    if anchor not in before:
        return {}
    ax, ay = before[anchor]
    rad = math.radians(delta_deg)
    cos_t = math.cos(rad)
    sin_t = math.sin(rad)
    rotated = {}
    for atom_id in side:
        if atom_id == anchor or atom_id not in before:
            continue
        x, y = before[atom_id]
        dx = x - ax
        dy = y - ay
        rotated[atom_id] = (ax + dx * cos_t - dy * sin_t, ay + dx * sin_t + dy * cos_t)
    return rotated


def _mirror_coords(
    coords: dict[int, tuple[float, float]],
    *,
    axis: str,
) -> dict[int, tuple[float, float]]:
    cx, cy = _center(coords)
    if axis == "y":
        return {atom_id: (2.0 * cx - x, y) for atom_id, (x, y) in coords.items()}
    return {atom_id: (x, 2.0 * cy - y) for atom_id, (x, y) in coords.items()}


def _has_explicit_stereo(graph: MolGraph, atom_ids: set[int], bonds: list[Bond]) -> bool:
    for atom_id in atom_ids:
        atom = graph.atoms.get(atom_id)
        if atom is None:
            continue
        if any(
            getattr(atom, attr, None)
            for attr in ("stereo_cip", "stereo_axial", "stereo_helical", "stereo_si_re")
        ):
            return True
    for bond in bonds:
        if any(
            getattr(bond, attr, None)
            for attr in ("stereo_ez", "stereo_axial", "stereo_endo_exo", "stereo_helical")
        ):
            return True
        stereo = getattr(getattr(bond, "stereo", None), "value", getattr(bond, "stereo", None))
        if stereo and str(stereo) not in {"none", "BondStereo.NONE"}:
            return True
    return False


def _ring_center_by_atom(
    ring_atoms: set[int],
    adjacency: dict[int, set[int]],
    coords: dict[int, tuple[float, float]],
) -> dict[int, tuple[float, float]]:
    centers: dict[int, tuple[float, float]] = {}
    core = _cycle_core_atoms(ring_atoms, adjacency)
    if not core:
        return centers
    center = _center({atom_id: coords[atom_id] for atom_id in core if atom_id in coords})
    for atom_id in core:
        centers[atom_id] = center
    return centers


def _tree_root(
    graph: MolGraph,
    component: set[int],
    adjacency: dict[int, set[int]],
    before: dict[int, tuple[float, float]],
) -> int:
    heavy = [
        atom_id
        for atom_id in component
        if graph.atoms.get(atom_id) is not None and graph.atoms[atom_id].element != "H"
    ]
    candidates = heavy or list(component)
    leaves = [atom_id for atom_id in candidates if len(adjacency.get(atom_id, set()) & component) <= 1]
    if leaves:
        candidates = leaves
    return min(candidates, key=lambda atom_id: (before.get(atom_id, (0.0, 0.0))[0], before.get(atom_id, (0.0, 0.0))[1], atom_id))


def _spread_angles(count: int, base: float) -> list[float]:
    if count <= 0:
        return []
    if count == 1:
        return [base]
    if count == 2:
        return [base + 30.0, base - 90.0]
    if count == 3:
        return [base, base + 120.0, base - 120.0]
    return [(base + idx * 360.0 / count) % 360.0 for idx in range(count)]


def _atom_geometry(graph: MolGraph, atom_id: int, bonds: list[Bond]) -> str:
    has_triple = False
    has_double = False
    for bond in bonds:
        if atom_id not in {bond.a1_id, bond.a2_id}:
            continue
        if int(getattr(bond, "order", 1) or 1) >= 3:
            has_triple = True
        if int(getattr(bond, "order", 1) or 1) == 2 or bool(getattr(bond, "is_aromatic", False)):
            has_double = True
    if has_triple:
        return "sp"
    if has_double:
        return "sp2"
    atom = graph.atoms.get(atom_id)
    if atom is not None and atom.element in {"B", "N", "O"} and len([b for b in bonds if atom_id in {b.a1_id, b.a2_id}]) <= 2:
        return "sp2"
    return "sp3"


def _bond_between(bonds: list[Bond], a1: int, a2: int) -> Bond | None:
    pair = {a1, a2}
    for bond in bonds:
        if {bond.a1_id, bond.a2_id} == pair:
            return bond
    return None


def _target_length_for_bond(bond: Bond | None, target: float) -> float:
    if bond is None:
        return target
    order = int(getattr(bond, "order", 1) or 1)
    if order >= 3:
        return target * 0.92
    if order == 2 or bool(getattr(bond, "is_aromatic", False)):
        return target * 0.97
    return target


def _rescale_to_bond_length(
    coords: dict[int, tuple[float, float]],
    bonds: list[Bond],
    target: float,
) -> dict[int, tuple[float, float]]:
    if not coords:
        return {}
    current = _average_bond_length(coords, bonds)
    if current <= 1e-6:
        return dict(coords)
    scale = target / current
    if not math.isfinite(scale):
        return dict(coords)
    scale = max(0.05, min(200.0, scale))
    cx, cy = _center(coords)
    return {
        atom_id: (cx + (x - cx) * scale, cy + (y - cy) * scale)
        for atom_id, (x, y) in coords.items()
    }


def _align_to_reference(
    reference: dict[int, tuple[float, float]],
    coords: dict[int, tuple[float, float]],
) -> dict[int, tuple[float, float]]:
    common = [atom_id for atom_id in reference if atom_id in coords]
    if not common:
        return dict(coords)
    ref_c = _center({atom_id: reference[atom_id] for atom_id in common})
    src_c = _center({atom_id: coords[atom_id] for atom_id in common})
    if len(common) == 1:
        return {
            atom_id: (x + ref_c[0] - src_c[0], y + ref_c[1] - src_c[1])
            for atom_id, (x, y) in coords.items()
        }

    def candidate(mirror_x: bool) -> tuple[dict[int, tuple[float, float]], float]:
        sum_cos = 0.0
        sum_sin = 0.0
        for atom_id in common:
            qx = coords[atom_id][0] - src_c[0]
            qy = coords[atom_id][1] - src_c[1]
            px = reference[atom_id][0] - ref_c[0]
            py = reference[atom_id][1] - ref_c[1]
            if mirror_x:
                qx = -qx
            sum_cos += qx * px + qy * py
            sum_sin += qx * py - qy * px
        theta = math.atan2(sum_sin, sum_cos)
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)
        transformed = {}
        error = 0.0
        for atom_id, (x, y) in coords.items():
            qx = x - src_c[0]
            qy = y - src_c[1]
            if mirror_x:
                qx = -qx
            rx = qx * cos_t - qy * sin_t + ref_c[0]
            ry = qx * sin_t + qy * cos_t + ref_c[1]
            transformed[atom_id] = (rx, ry)
            if atom_id in reference:
                error += (rx - reference[atom_id][0]) ** 2 + (ry - reference[atom_id][1]) ** 2
        return transformed, error

    direct, direct_error = candidate(False)
    mirrored, mirrored_error = candidate(True)
    return mirrored if mirrored_error < direct_error else direct


def _translate_to_center(
    coords: dict[int, tuple[float, float]],
    target_center: tuple[float, float],
) -> dict[int, tuple[float, float]]:
    if not coords:
        return {}
    cx, cy = _center(coords)
    return {atom_id: (x + target_center[0] - cx, y + target_center[1] - cy) for atom_id, (x, y) in coords.items()}


def _center(coords: dict[Any, tuple[float, float]]) -> tuple[float, float]:
    if not coords:
        return 0.0, 0.0
    return (
        sum(float(x) for x, _y in coords.values()) / len(coords),
        sum(float(y) for _x, y in coords.values()) / len(coords),
    )


def _average_bond_length(coords: dict[int, tuple[float, float]], bonds: list[Bond]) -> float:
    lengths = [
        _distance(coords[bond.a1_id], coords[bond.a2_id])
        for bond in bonds
        if bond.a1_id in coords and bond.a2_id in coords
    ]
    return sum(lengths) / len(lengths) if lengths else 0.0


def _bond_stats(coords: dict[int, tuple[float, float]], bonds: list[Bond]) -> dict[str, float]:
    lengths = [
        _distance(coords[bond.a1_id], coords[bond.a2_id])
        for bond in bonds
        if bond.a1_id in coords and bond.a2_id in coords
    ]
    if not lengths:
        return {"mean": 0.0, "min": 0.0, "max": 0.0}
    return {"mean": sum(lengths) / len(lengths), "min": min(lengths), "max": max(lengths)}


def _mean_displacement(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> float:
    distances = [
        _distance(before[atom_id], after[atom_id])
        for atom_id in atom_ids
        if atom_id in before and atom_id in after
    ]
    return sum(distances) / len(distances) if distances else 0.0


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(float(a[0]) - float(b[0]), float(a[1]) - float(b[1]))


def _endpoint(origin: tuple[float, float], angle_deg: float, length: float) -> tuple[float, float]:
    rad = math.radians(angle_deg)
    return (origin[0] + math.cos(rad) * length, origin[1] + math.sin(rad) * length)


def _angle_deg(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 360.0


def _angle_between(
    center: tuple[float, float],
    left: tuple[float, float],
    right: tuple[float, float],
) -> float:
    a1 = _angle_deg(center, left)
    a2 = _angle_deg(center, right)
    return _angle_distance(a1, a2)


def _angle_distance(a1: float, a2: float) -> float:
    return abs((a2 - a1 + 180.0) % 360.0 - 180.0)


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


def _bbox_ratio(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> float:
    def diag(coords: dict[int, tuple[float, float]]) -> float:
        xs = [coords[atom_id][0] for atom_id in atom_ids if atom_id in coords]
        ys = [coords[atom_id][1] for atom_id in atom_ids if atom_id in coords]
        if not xs or not ys:
            return 1.0
        return math.hypot((max(xs) - min(xs)) or 1.0, (max(ys) - min(ys)) or 1.0)

    before_diag = diag(before)
    after_diag = diag(after)
    if before_diag <= 1e-6:
        return 1.0
    return after_diag / before_diag
