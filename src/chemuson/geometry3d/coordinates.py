from __future__ import annotations

"""Utilidades de coordenadas 3D independientes de backends químicos."""

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def planar_coordinate_set(graph: MolGraph, *, source: str = "2d-fallback") -> CoordinateSet3D:
    """Crea una conformación plana usando las coordenadas actuales del canvas."""
    return CoordinateSet3D(
        positions={int(atom.id): (float(atom.x), float(atom.y), 0.0) for atom in graph.atoms.values()},
        source=source,
        method="planar-2d",
        metadata={"fallback": True},
    )


def coordinate_set_from_positions(
    positions: dict[int, tuple[float, float, float]],
    *,
    source: str,
    method: str,
    energy: float | None = None,
    metadata: dict[str, object] | None = None,
) -> CoordinateSet3D:
    """Construye un ``CoordinateSet3D`` normalizando IDs y floats."""
    return CoordinateSet3D(
        positions={int(k): (float(v[0]), float(v[1]), float(v[2])) for k, v in positions.items()},
        source=str(source),
        method=str(method),
        energy=energy,
        metadata=dict(metadata or {}),
    )
