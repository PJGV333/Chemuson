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

    @property
    def ok(self) -> bool:
        return self.coordinates is not None and bool(self.coordinates.positions)
