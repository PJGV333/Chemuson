from __future__ import annotations

"""Utilidades compartidas para exportadores computacionales."""

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def geometry_lines(graph: MolGraph, coordinates: CoordinateSet3D) -> list[str]:
    """Devuelve líneas `Elemento x y z` en orden de átomo."""
    positions = coordinates.normalized_positions()
    lines: list[str] = []
    for atom in sorted(graph.atoms.values(), key=lambda item: item.id):
        coords = positions.get(int(atom.id))
        if coords is None:
            continue
        x, y, z = coords
        lines.append(f"{atom.element:<3} {x: .8f} {y: .8f} {z: .8f}")
    return lines


def calculation_keyword(calculation: str) -> str:
    text = str(calculation or "sp").strip().lower()
    if text in {"opt", "optimization", "optimise", "optimize"}:
        return "opt"
    if text in {"freq", "frequency", "frequencies"}:
        return "freq"
    if text in {"optfreq", "opt+freq", "opt_freq"}:
        return "opt freq"
    return "sp"
