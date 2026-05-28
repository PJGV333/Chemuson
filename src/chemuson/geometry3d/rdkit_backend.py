from __future__ import annotations

"""Backend RDKit aislado para conformaciones y optimización 3D."""

from collections.abc import Iterator

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D, OptimizationFrame, OptimizationResult, OptimizationSettings


def generate_conformer(graph: MolGraph, settings: OptimizationSettings | None = None) -> OptimizationResult:
    """Genera una conformación 3D con RDKit en subproceso aislado."""
    settings = settings or OptimizationSettings()
    try:
        from chemuson.chemio.rdkit_safe import conformer_3d_isolated

        positions, metadata, error = conformer_3d_isolated(graph, timeout_s=settings.timeout_s)
    except Exception as exc:
        return OptimizationResult(None, method="rdkit", message=str(exc))
    if error or not positions:
        return OptimizationResult(None, method="rdkit", message=str(error or "empty_conformer"))
    coordset = CoordinateSet3D(
        positions=positions,
        source=str(metadata.get("source", "rdkit")),
        method=str(metadata.get("method", "rdkit")),
        energy=_float_or_none(metadata.get("energy")),
        metadata=metadata,
    )
    return OptimizationResult(coordset, converged=True, energy=coordset.energy, method=coordset.method)


def optimize(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
) -> OptimizationResult:
    """Optimiza geometría con UFF/MMFF usando RDKit aislado."""
    settings = settings or OptimizationSettings()
    try:
        from chemuson.chemio.rdkit_safe import optimize_3d_isolated

        positions, metadata, error = optimize_3d_isolated(
            graph,
            forcefield=settings.forcefield.value,
            max_iters=settings.max_iters,
            steps_per_update=settings.steps_per_update,
            timeout_s=settings.timeout_s,
        )
    except Exception as exc:
        return OptimizationResult(None, method="rdkit", message=str(exc))
    if error or not positions:
        return OptimizationResult(None, method="rdkit", message=str(error or "empty_optimization"))
    coordset = CoordinateSet3D(
        positions=positions,
        source=str(metadata.get("source", "rdkit")),
        method=str(metadata.get("method", settings.forcefield.value)),
        energy=_float_or_none(metadata.get("energy")),
        metadata=metadata,
    )
    converged = bool(metadata.get("converged", False))
    return OptimizationResult(
        coordset,
        converged=converged,
        energy=coordset.energy,
        method=coordset.method,
        message=str(metadata.get("message", "")),
    )


def optimize_iter(
    graph: MolGraph,
    coordset: CoordinateSet3D | None = None,
    settings: OptimizationSettings | None = None,
) -> Iterator[OptimizationFrame]:
    """Interfaz iterativa preparada para UI; RDKit actual devuelve frame final."""
    result = optimize(graph, coordset, settings)
    if result.coordinates is not None:
        yield OptimizationFrame(
            step=(settings or OptimizationSettings()).max_iters,
            coordinates=result.coordinates,
            energy=result.energy,
            converged=result.converged,
            message=result.message,
        )


def _float_or_none(value: object) -> float | None:
    try:
        return None if value is None else float(value)
    except Exception:
        return None
