from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable

from chemuson.chemio.depiction_candidates import block_donut_score, score_imported_depiction
from chemuson.clean2d.block_unwrap import block_unwrap_layout
from chemuson.clean2d.geometry import finite_coords, graph_atom_coords, graph_with_coords, normalize_atom_ids
from chemuson.clean2d.local_graph_cleaner import stereo_layout_signature
from chemuson.core.layers import BlockKind, build_multilayer_chemical_graph
from chemuson.core.model import MolGraph


@dataclass(frozen=True)
class ScaffoldDepictionReport:
    ok: bool
    reason: str
    source: str = "chemuson_scaffold_depiction"
    candidate_count: int = 0
    selected_strategy: str = ""
    score_before: float = 0.0
    score_after: float = 0.0
    donut_score_before: float = 0.0
    donut_score_after: float = 0.0
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ScaffoldLayoutCandidate:
    strategy: str
    coords: dict[int, tuple[float, float]]
    score: float
    donut_score: float
    rejected: bool = False
    rejection_reason: str = ""
    metadata: dict[str, object] = field(default_factory=dict)


def scaffold_depiction_candidates(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    target_bond_length: float = 40.0,
) -> list[ScaffoldLayoutCandidate]:
    selected = normalize_atom_ids(graph, atom_ids)
    if not selected:
        return []
    if not _is_complex_enough(graph, selected, target_bond_length):
        return [ScaffoldLayoutCandidate("not_complex", {}, math.inf, 0.0, True, "not_complex")]
    base_coords = graph_atom_coords(graph, selected)
    base_score, _ = score_imported_depiction(graph_with_coords(graph, base_coords), target_bond_length=target_bond_length)
    base_donut, _ = block_donut_score(graph_with_coords(graph, base_coords), target_bond_length=target_bond_length)

    raw: list[tuple[str, dict[int, tuple[float, float]] | None, dict[str, object]]] = []
    unwrap_coords, unwrap_report = block_unwrap_layout(graph, selected, target_bond_length=target_bond_length)
    raw.append(("longest_path_zigzag", unwrap_coords, {"unwrap_report": dict(unwrap_report.__dict__)}))
    if unwrap_coords:
        raw.extend(
            [
                ("preserve_rdkit_blocks", _affine_variant(unwrap_coords, x_scale=1.0, y_scale=1.18, shear=0.06), {}),
                ("force_block_layout", _affine_variant(unwrap_coords, x_scale=1.12, y_scale=0.92, shear=-0.10), {}),
                ("layered_scaffold", _affine_variant(unwrap_coords, x_scale=1.22, y_scale=1.08, shear=0.18), {}),
                ("radial_break", _affine_variant(unwrap_coords, x_scale=1.35, y_scale=0.82, shear=-0.16), {}),
            ]
        )
    candidates: list[ScaffoldLayoutCandidate] = []
    for strategy, coords, metadata in raw:
        if coords is None:
            candidates.append(ScaffoldLayoutCandidate(strategy, {}, math.inf, base_donut, True, "strategy_unavailable", metadata))
            continue
        rejection = _rejection_reason(graph, selected, base_coords, coords, target_bond_length, base_score, base_donut)
        candidate_graph = graph_with_coords(graph, coords)
        score, score_meta = score_imported_depiction(candidate_graph, target_bond_length=target_bond_length)
        donut, donut_meta = block_donut_score(candidate_graph, target_bond_length=target_bond_length)
        candidates.append(
            ScaffoldLayoutCandidate(
                strategy=strategy,
                coords=coords,
                score=score,
                donut_score=donut,
                rejected=bool(rejection),
                rejection_reason=rejection,
                metadata={**metadata, **score_meta, **donut_meta, "score_before": base_score, "donut_score_before": base_donut},
            )
        )
    candidates.sort(key=lambda item: (item.rejected, item.score + item.donut_score * 500.0, item.strategy))
    return candidates


def scaffold_depiction_layout(
    graph: MolGraph,
    atom_ids: Iterable[int] | None = None,
    *,
    target_bond_length: float = 40.0,
) -> tuple[dict[int, tuple[float, float]] | None, ScaffoldDepictionReport]:
    selected = normalize_atom_ids(graph, atom_ids)
    before = graph_atom_coords(graph, selected)
    score_before, _ = score_imported_depiction(graph_with_coords(graph, before), target_bond_length=target_bond_length) if before else (0.0, {})
    donut_before, _ = block_donut_score(graph_with_coords(graph, before), target_bond_length=target_bond_length) if before else (0.0, {})
    candidates = scaffold_depiction_candidates(graph, selected, target_bond_length=target_bond_length)
    accepted = [candidate for candidate in candidates if not candidate.rejected and candidate.coords]
    if not accepted:
        return None, ScaffoldDepictionReport(False, "no_safe_scaffold_candidate", candidate_count=len(candidates), score_before=score_before, donut_score_before=donut_before, metadata={"candidates": [_candidate_metadata(candidate) for candidate in candidates]})
    best = accepted[0]
    return best.coords, ScaffoldDepictionReport(True, "ok", candidate_count=len(candidates), selected_strategy=best.strategy, score_before=score_before, score_after=best.score, donut_score_before=donut_before, donut_score_after=best.donut_score, metadata={"candidates": [_candidate_metadata(candidate) for candidate in candidates]})


def _is_complex_enough(graph: MolGraph, selected: set[int], target: float) -> bool:
    try:
        layers = build_multilayer_chemical_graph(graph, selected)
    except Exception:
        return False
    counts: dict[BlockKind, int] = {}
    for block in layers.block_graph.blocks:
        counts[block.kind] = counts.get(block.kind, 0) + 1
    ring_count = sum(1 for motif in layers.motif_graph.motifs if getattr(motif, "kind", None).value == "ring")
    donut, _ = block_donut_score(graph, target_bond_length=target)
    return (
        counts.get(BlockKind.MACROCYCLE, 0) > 0
        or counts.get(BlockKind.CYCLOPHANE, 0) > 0
        or counts.get(BlockKind.AROMATIC_RING, 0) >= 4
        or ring_count >= 4
        or (len(selected) >= 35 and ring_count >= 3)
        or donut >= 4.0
    )


def _rejection_reason(graph: MolGraph, selected: set[int], before: dict[int, tuple[float, float]], after: dict[int, tuple[float, float]], target: float, score_before: float, donut_before: float) -> str:
    if set(after) != selected:
        return "missing_coordinates"
    if not finite_coords(after):
        return "non_finite_coordinates"
    after_graph = graph_with_coords(graph, after)
    score_after, _ = score_imported_depiction(after_graph, target_bond_length=target)
    donut_after, _ = block_donut_score(after_graph, target_bond_length=target)
    if not (score_after < score_before or donut_after < donut_before * 0.90):
        return "score_not_improved"
    try:
        if _wedge_count(graph) and stereo_layout_signature(graph, before, selected) != stereo_layout_signature(graph, after, selected):
            return "stereo_signature_changed"
    except Exception:
        pass
    return ""


def _wedge_count(graph: MolGraph) -> int:
    return sum(1 for bond in graph.bonds.values() if str(getattr(getattr(bond, "style", None), "value", getattr(bond, "style", ""))) in {"wedge", "hashed"})


def _affine_variant(coords: dict[int, tuple[float, float]], *, x_scale: float, y_scale: float, shear: float) -> dict[int, tuple[float, float]]:
    cx = sum(x for x, _y in coords.values()) / len(coords)
    cy = sum(y for _x, y in coords.values()) / len(coords)
    out: dict[int, tuple[float, float]] = {}
    for atom_id, (x, y) in coords.items():
        dx, dy = x - cx, y - cy
        out[atom_id] = (cx + dx * x_scale + dy * shear, cy + dy * y_scale)
    return out


def _candidate_metadata(candidate: ScaffoldLayoutCandidate) -> dict[str, object]:
    return {
        "strategy": candidate.strategy,
        "score": candidate.score,
        "donut_score": candidate.donut_score,
        "rejected": candidate.rejected,
        "rejection_reason": candidate.rejection_reason,
        "metadata": candidate.metadata,
    }
