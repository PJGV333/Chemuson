"""Explicit Haworth projection geometry for pyranose/furanose ring layouts."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

from chemuson.core.model import BondStyle

from .common import ProjectionGeometry, projection_atom, projection_bond

Orientation = Literal["up", "down"]


@dataclass(frozen=True, slots=True)
class HaworthSubstituent:
    """A declared group at a numbered ring carbon, projected up or down."""

    position: int
    label: str
    orientation: Orientation


@dataclass(frozen=True, slots=True)
class HaworthParameters:
    """Projection geometry and carbohydrate depiction policy, not renderer style."""

    ring_size: int = 6
    sugar_family: Literal["glucose", "generic"] = "glucose"
    bond_length: float = 40.0
    aspect_ratio: float = 1.0
    ring_skew: float = 0.0
    front_edge_policy: Literal["lower_arc", "none"] = "lower_arc"
    front_edge_style: Literal["bold", "plain"] = "bold"
    anomer: Literal["alpha", "beta"] = "beta"
    configuration: Literal["D", "L"] = "D"
    anomeric_up_override: bool | None = None
    substituents: tuple[HaworthSubstituent, ...] | None = None
    substituent_bond_length_ratio: float = 1.15

    def __post_init__(self) -> None:
        if self.ring_size not in (5, 6):
            raise ValueError("Haworth ring_size must be 5 (furanose) or 6 (pyranose)")
        if self.sugar_family not in ("glucose", "generic"):
            raise ValueError("sugar_family must be 'glucose' or 'generic'")
        for name, value in (
            ("bond_length", self.bond_length),
            ("aspect_ratio", self.aspect_ratio),
            ("ring_skew", self.ring_skew),
            ("substituent_bond_length_ratio", self.substituent_bond_length_ratio),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.bond_length <= 0.0:
            raise ValueError("bond_length must be positive")
        if not 0.5 <= self.aspect_ratio <= 2.0:
            raise ValueError("aspect_ratio must be between 0.5 and 2.0")
        if not -0.5 <= self.ring_skew <= 0.5:
            raise ValueError("ring_skew must be between -0.5 and 0.5")
        if self.front_edge_policy not in ("lower_arc", "none"):
            raise ValueError("front_edge_policy must be 'lower_arc' or 'none'")
        if self.front_edge_style not in ("bold", "plain"):
            raise ValueError("front_edge_style must be 'bold' or 'plain'")
        if self.anomer not in ("alpha", "beta"):
            raise ValueError("anomer must be 'alpha' or 'beta'")
        if self.configuration not in ("D", "L"):
            raise ValueError("configuration must be 'D' or 'L'")
        if self.anomeric_up_override is not None and not isinstance(
            self.anomeric_up_override, bool
        ):
            raise ValueError("anomeric_up_override must be a boolean or None")
        if self.substituent_bond_length_ratio <= 0.5:
            raise ValueError("substituent_bond_length_ratio must exceed 0.5")
        if self.substituents is not None:
            positions = [substituent.position for substituent in self.substituents]
            if len(positions) != len(set(positions)):
                raise ValueError("Haworth substituent positions must be unique")
            for substituent in self.substituents:
                if not 1 <= substituent.position < self.ring_size:
                    raise ValueError("substituent position must name a ring carbon")
                if not substituent.label.strip():
                    raise ValueError("substituent labels must not be empty")
                if substituent.orientation not in ("up", "down"):
                    raise ValueError("substituent orientation must be 'up' or 'down'")


def _glucose_substituents(parameters: HaworthParameters) -> tuple[HaworthSubstituent, ...]:
    """Return the declared D/L glucose pattern without assigning bond stereo."""
    is_d = parameters.configuration == "D"
    anomeric_up = parameters.anomeric_up_override
    if anomeric_up is None:
        anomeric_up = (parameters.anomer == "beta") == is_d
    down: Orientation = "down" if is_d else "up"
    up: Orientation = "up" if is_d else "down"
    return (
        HaworthSubstituent(1, "OH", "up" if anomeric_up else "down"),
        HaworthSubstituent(2, "OH", down),
        HaworthSubstituent(3, "OH", up),
        HaworthSubstituent(4, "OH", down),
        HaworthSubstituent(5, "CH2OH", up),
    )


def _ring_points(parameters: HaworthParameters) -> tuple[tuple[float, float], ...]:
    count = parameters.ring_size
    # The six-membered default is the conventional O-upper-right Haworth view.
    # A five-membered ring is similarly rotated to keep O on its upper-right arc.
    start_degrees = -60.0 if count == 6 else -54.0
    radius = parameters.bond_length / (2.0 * math.sin(math.pi / count))
    points = []
    for index in range(count):
        angle = math.radians(start_degrees + index * 360.0 / count)
        raw_x = radius * math.cos(angle)
        raw_y = radius * math.sin(angle)
        points.append(
            (
                raw_x * parameters.aspect_ratio + parameters.ring_skew * raw_y,
                raw_y,
            )
        )
    return tuple(points)


def _front_edges(
    points: tuple[tuple[float, float], ...], policy: str
) -> frozenset[int]:
    if policy == "none":
        return frozenset()
    center_y = sum(y for _x, y in points) / len(points)
    height = max(y for _x, y in points) - min(y for _x, y in points)
    threshold = center_y + height * 0.02
    return frozenset(
        index
        for index, point in enumerate(points)
        if (point[1] + points[(index + 1) % len(points)][1]) / 2.0 > threshold
    )


def _substituent_direction(
    point: tuple[float, float], orientation: Orientation
) -> tuple[float, float]:
    """Choose a stable diagonal, biasing inward when up/down opposes the arc."""
    length = math.hypot(*point)
    radial_x = point[0] / length
    radial_y = point[1] / length
    vertical = -1.0 if orientation == "up" else 1.0
    horizontal_sign = 1.0 if radial_x >= 0.0 else -1.0
    if radial_y * vertical < 0.0:
        horizontal_sign *= -1.0
    dx = horizontal_sign * abs(radial_x) * 0.42
    dy = vertical * math.sqrt(max(0.0, 1.0 - dx * dx))
    return dx, dy


def build_haworth_geometry(parameters: HaworthParameters | None = None) -> ProjectionGeometry:
    """Build explicit ring/front-arc/substituent geometry deterministically."""
    config = parameters or HaworthParameters()
    points = _ring_points(config)
    atom_specs = [projection_atom("O", *points[0], explicit=True)]
    atom_specs.extend(projection_atom("C", x, y, explicit=False) for x, y in points[1:])
    bond_specs = []
    front = _front_edges(points, config.front_edge_policy)
    bold = config.front_edge_style == "bold"
    for index in range(config.ring_size):
        edge_style = BondStyle.BOLD if bold and index in front else BondStyle.PLAIN
        bond_specs.append(projection_bond(index, (index + 1) % config.ring_size, style=edge_style))

    if config.substituents is not None:
        substituents = config.substituents
    elif config.sugar_family == "glucose" and config.ring_size == 6:
        substituents = _glucose_substituents(config)
    else:
        # Generic rings carry no inferred substituent/stereo information.
        substituents = ()

    for substituent in substituents:
        ring_index = substituent.position
        x, y = points[ring_index]
        dx, dy = _substituent_direction((x, y), substituent.orientation)
        if substituent.label == "CH2OH":
            dx = 0.0
            dy = -1.0 if substituent.orientation == "up" else 1.0
        endpoint_x = x + dx * config.bond_length * config.substituent_bond_length_ratio
        endpoint_y = y + dy * config.bond_length * config.substituent_bond_length_ratio
        substituent_index = len(atom_specs)
        atom_specs.append(projection_atom(substituent.label, endpoint_x, endpoint_y, explicit=True))
        bond_specs.append(projection_bond(ring_index, substituent_index))

    return ProjectionGeometry(
        family="haworth",
        bond_length=config.bond_length,
        atoms=tuple(atom_specs),
        bonds=tuple(bond_specs),
    )
