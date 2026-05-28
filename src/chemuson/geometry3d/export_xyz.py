from __future__ import annotations

"""Exportación XYZ simple para conformaciones 3D."""

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D


def molgraph_to_xyz(graph: MolGraph, coordinates: CoordinateSet3D, *, comment: str | None = None) -> str:
    """Exporta átomos con coordenadas 3D a formato XYZ."""
    positions = coordinates.normalized_positions()
    atoms = [atom for atom in sorted(graph.atoms.values(), key=lambda item: item.id) if int(atom.id) in positions]
    lines = [str(len(atoms)), comment if comment is not None else f"ChemUSON {coordinates.source} {coordinates.method}"]
    for atom in atoms:
        x, y, z = positions[int(atom.id)]
        lines.append(f"{atom.element} {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines) + "\n"


def xyz_to_positions(xyz: str, atom_ids: list[int] | tuple[int, ...]) -> dict[int, tuple[float, float, float]]:
    """Lee posiciones XYZ y las asigna por orden a IDs de átomo existentes."""
    lines = [line.strip() for line in str(xyz or "").splitlines() if line.strip()]
    if len(lines) < 2:
        return {}
    try:
        count = int(lines[0])
    except Exception:
        return {}
    data = lines[2 : 2 + count]
    positions: dict[int, tuple[float, float, float]] = {}
    for atom_id, line in zip(atom_ids, data):
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            positions[int(atom_id)] = (float(parts[1]), float(parts[2]), float(parts[3]))
        except Exception:
            continue
    return positions
