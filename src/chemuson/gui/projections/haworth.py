"""Explicit Haworth geometry and semantics for pyranose/furanose projections."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

from chemuson.core.model import BondStyle

from .common import ProjectionGeometry, projection_atom, projection_bond
from .semantics import HaworthProjectionSemantics, HaworthSubstituentSemantics

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
    perspective: float = 0.06
    tilt_degrees: float = 0.0
    front_edge_policy: Literal["lower_arc", "none"] = "lower_arc"
    front_edge_indices: tuple[int, ...] | None = None
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
            ("perspective", self.perspective),
            ("tilt_degrees", self.tilt_degrees),
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
        if not -0.25 <= self.perspective <= 0.25:
            raise ValueError("perspective must be between -0.25 and 0.25")
        if not -35.0 <= self.tilt_degrees <= 35.0:
            raise ValueError("tilt_degrees must be between -35 and 35")
        if self.front_edge_policy not in ("lower_arc", "none"):
            raise ValueError("front_edge_policy must be 'lower_arc' or 'none'")
        if self.front_edge_style not in ("bold", "plain"):
            raise ValueError("front_edge_style must be 'bold' or 'plain'")
        if self.front_edge_indices is not None:
            if any(
                not isinstance(index, int) or isinstance(index, bool)
                or not 0 <= index < self.ring_size
                for index in self.front_edge_indices
            ):
                raise ValueError("front_edge_indices must name ring edges")
            if len(set(self.front_edge_indices)) != len(self.front_edge_indices):
                raise ValueError("front_edge_indices must be unique")
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


# Vertices are numbered O, C1..C4 (furanose) or O, C1..C5 (pyranose).
# These deliberately non-regular outlines use near-uniform edge lengths while
# keeping the conventional oxygen-at-upper-right reading direction.
_CANONICAL_RING_POINTS: dict[int, tuple[tuple[float, float], ...]] = {
    5: (
        (0.00, 0.00),
        (0.52, 0.85),
        (0.12, 1.80),
        (-0.84, 1.52),
        (-0.98, 0.50),
    ),
    6: (
        (0.00, 0.00),
        (0.643, 0.766),
        (0.143, 1.632),
        (-0.857, 1.632),
        (-1.500, 0.866),
        (-1.000, 0.000),
    ),
}

# Position-specific exterior directions. Vertical sign is the declared Haworth
# up/down relation; the horizontal component clears the ring and neighbor labels.
_SUBSTITUENT_DIRECTIONS: dict[int, dict[int, dict[Orientation, tuple[float, float]]]] = {
    5: {
        1: {"up": (0.70, -0.714), "down": (0.70, 0.714)},
        2: {"up": (0.45, -0.893), "down": (0.45, 0.893)},
        3: {"up": (-0.985, -0.174), "down": (-0.45, 0.893)},
        4: {"up": (-0.70, -0.714), "down": (-0.966, 0.259)},
    },
    6: {
        1: {"up": (0.70, -0.714), "down": (0.70, 0.714)},
        2: {"up": (0.45, -0.893), "down": (0.45, 0.893)},
        3: {"up": (-0.985, -0.174), "down": (-0.45, 0.893)},
        4: {"up": (-0.70, -0.714), "down": (-0.966, 0.259)},
        5: {"up": (0.0, -1.0), "down": (0.0, 1.0)}
    },
}


def _ring_points(parameters: HaworthParameters) -> tuple[tuple[float, float], ...]:
    raw = _CANONICAL_RING_POINTS[parameters.ring_size]
    center_x = sum(x for x, _y in raw) / len(raw)
    center_y = sum(y for _x, y in raw) / len(raw)
    edge_lengths = tuple(
        math.hypot(
            raw[(index + 1) % len(raw)][0] - point[0],
            raw[(index + 1) % len(raw)][1] - point[1],
        )
        for index, point in enumerate(raw)
    )
    scale = parameters.bond_length / (sum(edge_lengths) / len(edge_lengths))
    height = (max(y for _x, y in raw) - min(y for _x, y in raw)) or 1.0
    angle = math.radians(parameters.tilt_degrees)
    cosine, sine = math.cos(angle), math.sin(angle)
    points = []
    for raw_x, raw_y in raw:
        y = (raw_y - center_y) * scale
        normalized_y = (raw_y - center_y) / height
        x = (raw_x - center_x) * scale * parameters.aspect_ratio
        x *= 1.0 + parameters.perspective * normalized_y
        x += parameters.ring_skew * y
        points.append((x * cosine - y * sine, x * sine + y * cosine))
    return tuple(points)


def _front_edges(parameters: HaworthParameters) -> tuple[int, ...]:
    if parameters.front_edge_indices is not None:
        return parameters.front_edge_indices
    if parameters.front_edge_policy == "none":
        return ()
    return (1, 2, 3)


def _substituent_direction(position: int, orientation: Orientation, ring_size: int) -> tuple[float, float]:
    try:
        return _SUBSTITUENT_DIRECTIONS[ring_size][position][orientation]
    except KeyError as exc:
        raise ValueError(
            f"no Haworth exterior direction for ring position {position} in ring size {ring_size}"
        ) from exc


def build_haworth_geometry(parameters: HaworthParameters | None = None) -> ProjectionGeometry:
    """Build explicit ring/front-edge/substituent geometry and its semantics."""
    config = parameters or HaworthParameters()
    points = _ring_points(config)
    atom_specs = [projection_atom("O", *points[0], explicit=True)]
    atom_specs.extend(projection_atom("C", x, y, explicit=False) for x, y in points[1:])
    bond_specs = []
    front_edges = _front_edges(config)
    bold = config.front_edge_style == "bold"
    for index in range(config.ring_size):
        edge_style = BondStyle.BOLD if bold and index in front_edges else BondStyle.PLAIN
        bond_specs.append(
            projection_bond(index, (index + 1) % config.ring_size, style=edge_style)
        )

    if config.substituents is not None:
        substituents = config.substituents
    elif config.sugar_family == "glucose" and config.ring_size == 6:
        substituents = _glucose_substituents(config)
    else:
        # Generic rings carry no inferred substituent/stereo information.
        substituents = ()

    semantic_substituents = []
    for substituent in substituents:
        ring_index = substituent.position
        x, y = points[ring_index]
        dx, dy = _substituent_direction(
            substituent.position, substituent.orientation, config.ring_size
        )
        endpoint_x = x + dx * config.bond_length * config.substituent_bond_length_ratio
        endpoint_y = y + dy * config.bond_length * config.substituent_bond_length_ratio
        substituent_index = len(atom_specs)
        atom_specs.append(
            projection_atom(substituent.label, endpoint_x, endpoint_y, explicit=True)
        )
        bond_specs.append(projection_bond(ring_index, substituent_index))
        semantic_substituents.append(
            HaworthSubstituentSemantics(
                ring_position=substituent.position,
                ring_atom_index=ring_index,
                group_atom_index=substituent_index,
                label=substituent.label,
                orientation=substituent.orientation,
            )
        )

    is_glucopyranose = config.sugar_family == "glucose" and config.ring_size == 6
    semantics = HaworthProjectionSemantics(
        ring_size=config.ring_size,
        ring_atom_indices=tuple(range(config.ring_size)),
        anomer=config.anomer if is_glucopyranose else None,
        configuration=config.configuration if is_glucopyranose else None,
        anomeric_position=1 if is_glucopyranose else None,
        substituents=tuple(semantic_substituents),
        front_edge_indices=front_edges,
    )
    return ProjectionGeometry(
        family="haworth",
        bond_length=config.bond_length,
        atoms=tuple(atom_specs),
        bonds=tuple(bond_specs),
        semantics=semantics,
    )
