"""Extensible chair-conformation geometry with axial/equatorial vectors."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

from chemuson.core.model import BondStyle

from .common import ProjectionGeometry, projection_atom, projection_bond

ChairOrientation = Literal["up", "down"]
ChairPlacement = Literal["axial", "equatorial"]

DEFAULT_AXIAL_VECTORS: tuple[tuple[float, float], ...] = (
    (0.0, -1.0),
    (0.0, 1.0),
    (0.0, -1.0),
    (0.0, 1.0),
    (0.0, -1.0),
    (0.0, 1.0),
)
DEFAULT_EQUATORIAL_VECTORS: tuple[tuple[float, float], ...] = (
    (-0.8, 0.6),
    (-0.6, -0.8),
    (0.5, -0.86),
    (0.9, -0.4),
    (0.9, 0.4),
    (-0.8, 0.6),
)


@dataclass(frozen=True, slots=True)
class ChairSubstituent:
    """A group attached at a ring carbon using an axial/equatorial vector."""

    position: int
    label: str
    orientation: ChairOrientation = "up"
    placement: ChairPlacement = "axial"


@dataclass(frozen=True, slots=True)
class ChairParameters:
    """Chair conformation and projection-specific substituent placement rules."""

    ring_bond_length: float = 40.0
    conformation: Literal["chair_a", "chair_b"] = "chair_a"
    front_edge_policy: Literal["lower_arc", "none"] = "lower_arc"
    axial_vectors: tuple[tuple[float, float], ...] = DEFAULT_AXIAL_VECTORS
    equatorial_vectors: tuple[tuple[float, float], ...] = DEFAULT_EQUATORIAL_VECTORS
    axial_bond_length_ratio: float = 1.0
    equatorial_bond_length_ratio: float = 1.0
    substituents: tuple[ChairSubstituent, ...] = ()

    def __post_init__(self) -> None:
        if not math.isfinite(self.ring_bond_length) or self.ring_bond_length <= 0.0:
            raise ValueError("ring_bond_length must be finite and positive")
        if self.conformation not in ("chair_a", "chair_b"):
            raise ValueError("conformation must be 'chair_a' or 'chair_b'")
        if self.front_edge_policy not in ("lower_arc", "none"):
            raise ValueError("front_edge_policy must be 'lower_arc' or 'none'")
        for name, ratio in (
            ("axial_bond_length_ratio", self.axial_bond_length_ratio),
            ("equatorial_bond_length_ratio", self.equatorial_bond_length_ratio),
        ):
            if not math.isfinite(ratio) or ratio <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for vectors in (self.axial_vectors, self.equatorial_vectors):
            if len(vectors) != 6:
                raise ValueError("chair vector policies must contain six vectors")
            if any(
                not math.isfinite(x)
                or not math.isfinite(y)
                or math.hypot(x, y) == 0.0
                for x, y in vectors
            ):
                raise ValueError("chair vectors must be finite and non-zero")
        positions = [item.position for item in self.substituents]
        if len(positions) != len(set(positions)):
            raise ValueError("chair substituent positions must be unique")
        for item in self.substituents:
            if not 0 <= item.position < 6:
                raise ValueError("chair substituent position must be between 0 and 5")
            if not item.label.strip():
                raise ValueError("chair substituent labels must not be empty")
            if item.orientation not in ("up", "down"):
                raise ValueError("chair substituent orientation must be 'up' or 'down'")
            if item.placement not in ("axial", "equatorial"):
                raise ValueError("chair substituent placement must be axial or equatorial")


def _ring_coordinates(config: ChairParameters) -> tuple[tuple[float, float], ...]:
    base = (
        (-1.35, 0.35),
        (-0.80, -0.55),
        (0.15, -0.25),
        (1.10, -0.55),
        (0.55, 0.35),
        (-0.35, 0.05),
    )
    mirror = -1.0 if config.conformation == "chair_b" else 1.0
    return tuple((x * config.ring_bond_length, y * mirror * config.ring_bond_length) for x, y in base)


def _front_edges(points: tuple[tuple[float, float], ...], policy: str) -> frozenset[int]:
    if policy == "none":
        return frozenset()
    center_y = sum(y for _x, y in points) / len(points)
    height = max(y for _x, y in points) - min(y for _x, y in points)
    threshold = center_y + height * 0.05
    return frozenset(
        index
        for index, point in enumerate(points)
        if (point[1] + points[(index + 1) % len(points)][1]) / 2.0 > threshold
    )


def _oriented_vector(
    vector: tuple[float, float], orientation: ChairOrientation, mirror: float
) -> tuple[float, float]:
    x, y = vector[0], vector[1] * mirror
    if (orientation == "up" and y > 0.0) or (orientation == "down" and y < 0.0):
        x, y = -x, -y
    length = math.hypot(x, y)
    return x / length, y / length


def build_chair_geometry(parameters: ChairParameters | None = None) -> ProjectionGeometry:
    """Build chair ring and optional explicitly placed substituents."""
    config = parameters or ChairParameters()
    points = _ring_coordinates(config)
    atoms = [projection_atom("C", x, y, explicit=False) for x, y in points]
    bold_edges = _front_edges(points, config.front_edge_policy)
    bonds = [
        projection_bond(
            index,
            (index + 1) % 6,
            style=BondStyle.BOLD if index in bold_edges else BondStyle.PLAIN,
        )
        for index in range(6)
    ]
    mirror = -1.0 if config.conformation == "chair_b" else 1.0
    for substituent in config.substituents:
        vectors = (
            config.axial_vectors
            if substituent.placement == "axial"
            else config.equatorial_vectors
        )
        dx, dy = _oriented_vector(vectors[substituent.position], substituent.orientation, mirror)
        ratio = (
            config.axial_bond_length_ratio
            if substituent.placement == "axial"
            else config.equatorial_bond_length_ratio
        )
        atom_index = len(atoms)
        x, y = points[substituent.position]
        atoms.append(
            projection_atom(
                substituent.label,
                x + dx * config.ring_bond_length * ratio,
                y + dy * config.ring_bond_length * ratio,
                explicit=True,
            )
        )
        bonds.append(projection_bond(substituent.position, atom_index))
    return ProjectionGeometry(
        family="chair",
        bond_length=config.ring_bond_length,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
    )
