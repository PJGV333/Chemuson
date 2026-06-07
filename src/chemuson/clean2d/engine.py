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
from chemuson.clean2d.safety import (
    Clean2DQualityReport,
    evaluate_clean2d_layout,
    has_cycles,
    is_clean2d_candidate_safe,
    min_nonbonded_distance,
    ring_degeneracy_score,
)
from chemuson.clean2d.v2 import Clean2DParameters, optimize_clean2d_positions
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
    selected = _normalize_atom_ids(graph, atom_ids)
    before = _coords_for_atoms(graph, selected)
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
    candidates: list[Clean2DCandidate] = []

    current = _normalized_candidate_coords(before, before, bonds, target)
    candidates.append(
        Clean2DCandidate(
            source="current",
            coords=current,
            message="Estructura 2D ya estaba limpia",
        )
    )

    if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION}:
        candidates.extend(
            candidate
            for candidate in (
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
            "visual_score": candidate_quality.visual_score,
        }
        safety_mode = _safety_mode_for_candidate(mode, candidate, baseline_bad)
        safe = is_clean2d_candidate_safe(report, safety_mode)
        if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION} and candidate.source == "current":
            if quality.quality_class == "needs_rebuild":
                safe = False
                report.rejection_reason = "geometria_actual_necesita_reconstruccion"
            elif quality.quality_class == "needs_polish":
                safe = False
                report.rejection_reason = "geometria_actual_requiere_optimizacion"
        if mode in {Clean2DMode.QUICK, Clean2DMode.PUBLICATION} and candidate_quality.quality_class == "needs_rebuild":
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

        if candidate.source == "current" and quality.quality_class == "good":
            score = min(score, baseline_score * 0.25)
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

    if length_max > 0.14 or length_rms > 0.055:
        return _replace_quality(base_report, quality_class="needs_polish", reason="longitudes_suboptimas")
    if angle_stats["max"] > 34.0 or angle_stats["rms"] > 18.0 or angle_stats["severe_count"] > 0:
        return _replace_quality(base_report, quality_class="needs_polish", reason="angulos_suboptimos")
    if nonbonded < target * 0.50:
        return _replace_quality(base_report, quality_class="needs_polish", reason="distancias_no_enlazadas_estrechas")
    if min_ring != math.inf and min_ring < 0.25:
        return _replace_quality(base_report, quality_class="needs_polish", reason="anillo_suboptimo")
    return base_report


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
    coords = _normalized_candidate_coords(before, coords, bonds, target, align=False)
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
        coords = _place_ring_template(ring, before, out, placed, ring_centers, target)
        if not coords:
            continue
        for atom_id, coord in coords.items():
            out[atom_id] = coord
            placed.add(atom_id)
        ring_centers[ring_idx] = _center({atom_id: out[atom_id] for atom_id in ring})

    if placed:
        _layout_branches_from_placed_core(graph, atom_ids, bonds, adjacency, out, placed, target)
        for ring in pending_rings:
            if not (set(ring) & placed):
                continue
            coords = _place_ring_template(ring, before, out, placed, ring_centers, target)
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
) -> dict[int, tuple[float, float]]:
    shared_edges: list[tuple[int, int, int]] = []
    n = len(ring)
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
) -> None:
    ring_center_by_atom = _ring_center_by_atom(placed, adjacency, out)
    queue = sorted(placed)
    parent: dict[int, int | None] = {atom_id: None for atom_id in placed}
    while queue:
        anchor = queue.pop(0)
        unplaced = sorted(neigh for neigh in adjacency.get(anchor, set()) if neigh in atom_ids and neigh not in placed)
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
        "rdkit_isolated": 1,
        "rdkit_direct": 2,
        "internal_templates": 3,
        "clean2d_v2": 4,
        "safe_fallback": 5,
        "propose_reflection": 1,
        "propose_rotation": 2,
        "propose_3d_projection": 3,
        "propose_mirror": 4,
    }
    return priorities.get(source, 99)


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
