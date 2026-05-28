from __future__ import annotations

"""Modelos puros para conformaciones 3D y optimización molecular."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping


class ForceField(str, Enum):
    """Force fields soportados por backends externos."""

    UFF = "UFF"
    MMFF94 = "MMFF94"
    MMFF94S = "MMFF94s"
    GAFF = "GAFF"
    GHEMICAL = "Ghemical"


@dataclass(frozen=True)
class CoordinateSet3D:
    """Coordenadas 3D asociadas a IDs de átomos de un ``MolGraph``."""

    positions: Mapping[int, tuple[float, float, float]]
    source: str
    method: str
    energy: float | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def normalized_positions(self) -> dict[int, tuple[float, float, float]]:
        """Devuelve posiciones normalizadas a tipos primitivos."""
        return {
            int(atom_id): (float(coords[0]), float(coords[1]), float(coords[2]))
            for atom_id, coords in self.positions.items()
        }


@dataclass(frozen=True)
class OptimizationSettings:
    """Parámetros comunes para generación/optimización 3D."""

    forcefield: ForceField = ForceField.MMFF94
    max_iters: int = 200
    steps_per_update: int = 25
    timeout_s: float = 20.0
    convergence: float | None = None
    seed: int = 0xC0FFEE


@dataclass(frozen=True)
class OptimizationFrame:
    """Snapshot intermedio de una optimización molecular."""

    step: int
    coordinates: CoordinateSet3D
    energy: float | None = None
    converged: bool = False
    message: str = ""


@dataclass(frozen=True)
class OptimizationResult:
    """Resultado final de generación u optimización 3D."""

    coordinates: CoordinateSet3D | None
    converged: bool = False
    energy: float | None = None
    method: str = ""
    message: str = ""
    frames: tuple[OptimizationFrame, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.coordinates is not None and bool(self.coordinates.positions)


@dataclass(frozen=True)
class SceneAtom3D:
    """Átomo listo para una escena/visor 3D, incluyendo átomos generados."""

    id: int | str
    element: str
    x: float
    y: float
    z: float
    radius: float
    color: str
    generated: bool = False


@dataclass(frozen=True)
class SceneBond3D:
    """Enlace listo para render 3D."""

    id: int | str
    a1_id: int | str
    a2_id: int | str
    order: int = 1


@dataclass(frozen=True)
class SceneMolecule3D:
    """Representación desacoplada para un visor 3D futuro.

    ``MolGraph`` sigue siendo el modelo 2D editable. Esta escena puede incluir
    hidrógenos generados por backends externos sin mutar el documento 2D.
    """

    atoms: tuple[SceneAtom3D, ...]
    bonds: tuple[SceneBond3D, ...]
    source: str = ""


_VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}

_ELEMENT_COLORS = {
    "H": "#FFFFFF",
    "C": "#404040",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#90E050",
    "P": "#FF8000",
    "S": "#FFFF30",
    "Cl": "#1FF01F",
    "Br": "#A62929",
    "I": "#940094",
}


def scene_molecule_from_graph(graph, coordinates: CoordinateSet3D) -> SceneMolecule3D:
    """Construye una escena 3D desde ``MolGraph`` y coordenadas asociadas."""
    positions = coordinates.normalized_positions()
    atoms: list[SceneAtom3D] = []
    for atom in sorted(getattr(graph, "atoms", {}).values(), key=lambda item: item.id):
        pos = positions.get(int(atom.id))
        if pos is None:
            continue
        element = str(atom.element)
        atoms.append(
            SceneAtom3D(
                id=int(atom.id),
                element=element,
                x=pos[0],
                y=pos[1],
                z=pos[2],
                radius=float(_VDW_RADII.get(element, 1.70)),
                color=str(_ELEMENT_COLORS.get(element, "#B0B0B0")),
                generated=False,
            )
        )
    generated_atoms = coordinates.metadata.get("generated_atoms", ()) if hasattr(coordinates, "metadata") else ()
    if isinstance(generated_atoms, list):
        for item in generated_atoms:
            if not isinstance(item, Mapping):
                continue
            coords = item.get("coords") or item.get("position")
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                continue
            element = str(item.get("element", "H"))
            atoms.append(
                SceneAtom3D(
                    id=str(item.get("id", f"gen{len(atoms) + 1}")),
                    element=element,
                    x=float(coords[0]),
                    y=float(coords[1]),
                    z=float(coords[2]),
                    radius=float(_VDW_RADII.get(element, 1.20)),
                    color=str(_ELEMENT_COLORS.get(element, "#B0B0B0")),
                    generated=True,
                )
            )
    bonds = tuple(
        SceneBond3D(int(bond.id), int(bond.a1_id), int(bond.a2_id), int(getattr(bond, "order", 1) or 1))
        for bond in sorted(getattr(graph, "bonds", {}).values(), key=lambda item: item.id)
        if int(bond.a1_id) in positions and int(bond.a2_id) in positions
    )
    return SceneMolecule3D(atoms=tuple(atoms), bonds=bonds, source=coordinates.source)
