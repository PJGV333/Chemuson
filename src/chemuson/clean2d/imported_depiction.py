from __future__ import annotations

import copy
from dataclasses import dataclass, field
import json
import math
import os
from pathlib import Path
import sys
import tempfile

from chemuson.chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
    rdkit_worker_unavailable_message,
    smiles_to_molgraph_isolated_or_error,
)
from chemuson.chemio.rdkit_safe import smiles_depict_candidates_isolated
from chemuson.clean2d.block_unwrap import block_unwrap_layout
from chemuson.clean2d.depiction_quality import score_imported_depiction
from chemuson.clean2d.geometry import apply_coords_in_place
from chemuson.clean2d.scaffold_depiction import scaffold_depiction_candidates
from chemuson.core.model import MolGraph


@dataclass(frozen=True)
class DepictionCandidate:
    source: str
    graph: MolGraph
    score: float
    rejected: bool = False
    rejection_reason: str = ""
    metadata: dict[str, object] = field(default_factory=dict)


def smiles_to_depiction_candidates(
    smiles: str,
    *,
    target_bond_length: float = 40.0,
    timeout_s: float = 15.0,
) -> list[DepictionCandidate]:
    """Devuelve candidatos de depiction SMILES ordenados por calidad."""
    clean_smiles = str(smiles or "").strip()
    if not clean_smiles:
        raise ValueError("SMILES vacío")
    raw_candidates, error = smiles_depict_candidates_isolated(
        clean_smiles,
        target_bond_length=target_bond_length,
        timeout_s=timeout_s,
    )
    candidates: list[DepictionCandidate] = []
    for raw in raw_candidates:
        source = str(raw.get("source", "rdkit_candidate") or "rdkit_candidate")
        metadata = dict(raw.get("metadata", {}) if isinstance(raw.get("metadata"), dict) else {})
        if not raw.get("ok"):
            candidates.append(_rejected_depiction_candidate(source, str(raw.get("error", "candidate_failed") or "candidate_failed"), metadata))
            continue
        molblock = str(raw.get("molblock", "") or "")
        if not molblock.strip():
            candidates.append(_rejected_depiction_candidate(source, "empty_molblock", metadata))
            continue
        try:
            graph = molfile_to_molgraph(molblock, target_bond_length=target_bond_length)
            score, score_metadata = score_imported_depiction(graph, target_bond_length=target_bond_length)
            candidates.append(
                DepictionCandidate(
                    source=source,
                    graph=graph,
                    score=score,
                    metadata={**metadata, **score_metadata, "worker_error": error or ""},
                )
            )
            unwrap = _block_unwrap_depiction_candidate(graph, source, score, target_bond_length)
            if unwrap is not None:
                candidates.append(unwrap)
            candidates.extend(_scaffold_depiction_candidates(graph, source, score, target_bond_length))
        except Exception as exc:
            candidates.append(_rejected_depiction_candidate(source, str(exc), metadata))
    candidates.sort(key=lambda item: (item.rejected, item.score, item.source))
    if not candidates and error:
        raise RuntimeError(rdkit_worker_unavailable_message(str(error)))
    return candidates


def _rejected_depiction_candidate(source: str, reason: str, metadata: dict[str, object]) -> DepictionCandidate:
    return DepictionCandidate(source=source, graph=MolGraph(), score=math.inf, rejected=True, rejection_reason=reason, metadata=metadata)


def _scaffold_depiction_candidates(graph: MolGraph, parent_source: str, parent_score: float, target_bond_length: float) -> list[DepictionCandidate]:
    try:
        scaffold_candidates = scaffold_depiction_candidates(graph, target_bond_length=target_bond_length)
    except Exception:
        return []
    out: list[DepictionCandidate] = []
    for candidate in scaffold_candidates:
        if candidate.rejected or not candidate.coords:
            continue
        scaffold_graph = copy.deepcopy(graph)
        apply_coords_in_place(scaffold_graph, candidate.coords)
        score, metadata = score_imported_depiction(scaffold_graph, target_bond_length=target_bond_length)
        out.append(
            DepictionCandidate(
                source=f"chemuson_scaffold:{candidate.strategy}:{parent_source}",
                graph=scaffold_graph,
                score=score,
                metadata={
                    **metadata,
                    "parent_source": parent_source,
                    "parent_score": parent_score,
                    "strategy": candidate.strategy,
                    "scaffold_score": candidate.score,
                    "scaffold_donut_score": candidate.donut_score,
                    "stereo_wedge_count": _visual_wedge_count(scaffold_graph),
                    "scaffold_metadata": dict(candidate.metadata),
                },
            )
        )
    return out


def _block_unwrap_depiction_candidate(graph: MolGraph, parent_source: str, parent_score: float, target_bond_length: float) -> DepictionCandidate | None:
    try:
        coords, report = block_unwrap_layout(graph, target_bond_length=target_bond_length)
    except Exception:
        return None
    if coords is None or not report.ok:
        return None
    unwrapped = copy.deepcopy(graph)
    apply_coords_in_place(unwrapped, coords)
    unwrap_score, unwrap_metadata = score_imported_depiction(unwrapped, target_bond_length=target_bond_length)
    report_metadata = dict(report.__dict__)
    report_metadata.pop("metadata", None)
    return DepictionCandidate(
        source=f"chemuson_block_unwrap:{parent_source}",
        graph=unwrapped,
        score=unwrap_score,
        metadata={
            **unwrap_metadata,
            "unwrap_report": report_metadata,
            "unwrap_report_metadata": dict(report.metadata),
            "parent_source": parent_source,
            "parent_score": parent_score,
            "unwrap_score": unwrap_score,
            "donut_score_before": report.donut_score_before,
            "donut_score_after": report.donut_score_after,
        },
    )


def smiles_to_molgraph_best_depiction(smiles: str, *, target_bond_length: float = 40.0, timeout_s: float = 15.0) -> MolGraph:
    """Importa SMILES usando el mejor candidato de depiction aislado disponible."""
    graph, _report = smiles_to_molgraph_best_depiction_with_report(smiles, target_bond_length=target_bond_length, timeout_s=timeout_s)
    return graph


def smiles_to_molgraph_best_depiction_with_report(smiles: str, *, target_bond_length: float = 40.0, timeout_s: float = 15.0) -> tuple[MolGraph, dict[str, object]]:
    """Importa SMILES y devuelve también ranking/metadata de depiction."""
    clean_smiles = str(smiles or "").strip()
    if not clean_smiles:
        raise ValueError("SMILES vacío")
    candidate_error = ""
    try:
        candidates = smiles_to_depiction_candidates(clean_smiles, target_bond_length=target_bond_length, timeout_s=timeout_s)
        ranked = sorted(candidates, key=lambda item: (item.rejected, item.score, item.source))
        _debug_depiction_candidates(ranked)
        for candidate in ranked:
            if not candidate.rejected and candidate.graph.atoms:
                return candidate.graph, {"selected_source": candidate.source, "selected_score": candidate.score, "selected_metadata": dict(candidate.metadata), "candidates": [_depiction_candidate_row(item) for item in ranked]}
        rejected = [candidate.rejection_reason for candidate in candidates if candidate.rejected and candidate.rejection_reason]
        candidate_error = "; ".join(rejected) or "no_usable_depiction_candidates"
    except Exception as exc:
        candidate_error = str(exc)
    try:
        graph = smiles_to_molgraph_isolated_or_error(clean_smiles, timeout_s=min(timeout_s, 8.0))
        return graph, {"selected_source": "isolated_to_molblock_fallback", "error": candidate_error}
    except Exception as exc:
        fallback_error = str(exc)
    raise RuntimeError(rdkit_worker_unavailable_message(candidate_error or fallback_error, direct_error=fallback_error))


def _depiction_candidate_row(candidate: DepictionCandidate) -> dict[str, object]:
    return {"source": candidate.source, "score": candidate.score, "donut_score": candidate.metadata.get("donut_score", candidate.metadata.get("donut_score_after", "")), "rejected": candidate.rejected, "rejection_reason": candidate.rejection_reason, "atom_count": len(candidate.graph.atoms), "bond_count": len(candidate.graph.bonds), "wedge_count": _visual_wedge_count(candidate.graph), "metadata": dict(candidate.metadata)}


def _visual_wedge_count(graph: MolGraph) -> int:
    return sum(1 for bond in graph.bonds.values() if getattr(bond.style, "value", str(bond.style)) in {"wedge", "hashed"})


def _debug_depiction_candidates(candidates: list[DepictionCandidate]) -> None:
    if str(os.environ.get("CHEMUSON_DEBUG_DEPICTION", "")).strip().lower() not in {"1", "true", "yes", "on"}:
        return
    rows = [_depiction_candidate_row(candidate) for candidate in candidates]
    for idx, candidate in enumerate(candidates, start=1):
        status = "rejected" if candidate.rejected else "ok"
        donut = candidate.metadata.get("donut_score", candidate.metadata.get("donut_score_after", ""))
        sys.stderr.write(f"[chemuson depiction] #{idx} source={candidate.source} status={status} score={candidate.score:.3f} donut={donut} reason={candidate.rejection_reason}\n")
    try:
        folder = Path(tempfile.gettempdir()) / "chemuson_depiction_debug"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "ranking.json").write_text(json.dumps(rows, indent=2, default=str), encoding="utf-8")
        for idx, candidate in enumerate(candidates, start=1):
            try:
                (folder / f"candidate_{idx:02d}_{candidate.source.replace(':', '_')}.mol").write_text(molgraph_to_molfile(candidate.graph), encoding="utf-8")
            except Exception:
                pass
    except Exception:
        pass
