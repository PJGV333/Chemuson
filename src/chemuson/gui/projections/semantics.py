"""Typed semantic descriptors for graphical projection conventions.

These records describe projection intent in memory. They are deliberately
separate from MolGraph BondStereo and CMSN persistence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

ProjectionOrientation = Literal["up", "down"]
HaworthAnomer = Literal["alpha", "beta"]
CarbohydrateConfiguration = Literal["D", "L"]
FischerConfiguration = Literal["R", "S"]
ChairConformation = Literal["chair_a", "chair_b"]
ChairPlacement = Literal["axial", "equatorial"]


def _validate_indices(indices: tuple[int, ...], *, name: str) -> None:
    if any(not isinstance(index, int) or isinstance(index, bool) or index < 0 for index in indices):
        raise ValueError(f"{name} must contain non-negative atom indices")


@dataclass(frozen=True, slots=True)
class HaworthSubstituentSemantics:
    """A declared group and up/down relation at one numbered ring carbon."""

    ring_position: int
    ring_atom_index: int
    group_atom_index: int
    label: str
    orientation: ProjectionOrientation

    def __post_init__(self) -> None:
        _validate_indices((self.ring_atom_index, self.group_atom_index), name="Haworth indices")
        if self.ring_position < 1:
            raise ValueError("Haworth ring position must name a carbon")
        if not self.label.strip():
            raise ValueError("Haworth substituent label must not be empty")
        if self.orientation not in ("up", "down"):
            raise ValueError("Haworth orientation must be 'up' or 'down'")


@dataclass(frozen=True, slots=True)
class HaworthProjectionSemantics:
    """Alpha/beta, D/L and substituent semantics separate from coordinates."""

    ring_size: int
    ring_atom_indices: tuple[int, ...]
    anomer: HaworthAnomer | None
    configuration: CarbohydrateConfiguration | None
    anomeric_position: int | None
    substituents: tuple[HaworthSubstituentSemantics, ...]
    front_edge_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.ring_size not in (5, 6) or len(self.ring_atom_indices) != self.ring_size:
            raise ValueError("Haworth semantics require a five- or six-membered indexed ring")
        _validate_indices(self.ring_atom_indices, name="Haworth ring indices")
        if len(set(self.ring_atom_indices)) != self.ring_size:
            raise ValueError("Haworth ring atom indices must be unique")
        if self.anomer not in (None, "alpha", "beta"):
            raise ValueError("invalid Haworth anomer")
        if self.configuration not in (None, "D", "L"):
            raise ValueError("invalid carbohydrate configuration")
        if self.anomeric_position is not None and not 1 <= self.anomeric_position < self.ring_size:
            raise ValueError("invalid Haworth anomeric position")
        if any(not 0 <= edge < self.ring_size for edge in self.front_edge_indices):
            raise ValueError("Haworth front-edge index is outside the ring")
        if len(set(self.front_edge_indices)) != len(self.front_edge_indices):
            raise ValueError("Haworth front-edge indices must be unique")
        positions: set[int] = set()
        for substituent in self.substituents:
            if not 1 <= substituent.ring_position < self.ring_size:
                raise ValueError("Haworth substituent position is outside the ring")
            if substituent.ring_atom_index != self.ring_atom_indices[substituent.ring_position]:
                raise ValueError("Haworth substituent ring index does not match its position")
            if substituent.ring_position in positions:
                raise ValueError("Haworth semantic substituent positions must be unique")
            positions.add(substituent.ring_position)

    @property
    def referenced_atom_indices(self) -> tuple[int, ...]:
        return self.ring_atom_indices + tuple(
            index
            for substituent in self.substituents
            for index in (substituent.ring_atom_index, substituent.group_atom_index)
        )


@dataclass(frozen=True, slots=True)
class FischerCenterSemantics:
    """Indexed substituent order at one Fischer cross; R/S is never inferred."""

    center_atom_index: int
    left_atom_index: int
    right_atom_index: int
    upper_atom_index: int
    lower_atom_index: int
    absolute_configuration: FischerConfiguration | None = None

    def __post_init__(self) -> None:
        _validate_indices(
            (
                self.center_atom_index,
                self.left_atom_index,
                self.right_atom_index,
                self.upper_atom_index,
                self.lower_atom_index,
            ),
            name="Fischer indices",
        )
        if len(
            {
                self.center_atom_index,
                self.left_atom_index,
                self.right_atom_index,
                self.upper_atom_index,
                self.lower_atom_index,
            }
        ) != 5:
            raise ValueError("Fischer center and four substituent indices must be distinct")
        if self.absolute_configuration not in (None, "R", "S"):
            raise ValueError("Fischer absolute configuration must be R, S or undeclared")


@dataclass(frozen=True, slots=True)
class FischerProjectionSemantics:
    """The Fischer crossing convention and explicit indexed centers."""

    centers: tuple[FischerCenterSemantics, ...]
    horizontal_toward_observer: bool = True
    vertical_away_from_observer: bool = True

    def __post_init__(self) -> None:
        if not self.centers:
            raise ValueError("Fischer semantics require at least one stereocenter")
        if self.horizontal_toward_observer is not True:
            raise ValueError("Fischer horizontal bonds must point toward the observer")
        if self.vertical_away_from_observer is not True:
            raise ValueError("Fischer vertical bonds must recede from the observer")

    @property
    def referenced_atom_indices(self) -> tuple[int, ...]:
        return tuple(
            index
            for center in self.centers
            for index in (
                center.center_atom_index,
                center.left_atom_index,
                center.right_atom_index,
                center.upper_atom_index,
                center.lower_atom_index,
            )
        )


@dataclass(frozen=True, slots=True)
class ChairSubstituentSemantics:
    """One indexed substituent's declared up/down and axial/equatorial state."""

    ring_position: int
    ring_atom_index: int
    group_atom_index: int
    label: str
    orientation: ProjectionOrientation
    placement: ChairPlacement

    def __post_init__(self) -> None:
        _validate_indices((self.ring_atom_index, self.group_atom_index), name="chair indices")
        if not 0 <= self.ring_position < 6:
            raise ValueError("chair ring position must be between 0 and 5")
        if not self.label.strip():
            raise ValueError("chair substituent label must not be empty")
        if self.orientation not in ("up", "down"):
            raise ValueError("chair orientation must be 'up' or 'down'")
        if self.placement not in ("axial", "equatorial"):
            raise ValueError("chair placement must be axial or equatorial")


@dataclass(frozen=True, slots=True)
class ChairProjectionSemantics:
    """Conformation and explicit substituent relations for a six-membered chair."""

    conformation: ChairConformation
    ring_atom_indices: tuple[int, ...]
    substituents: tuple[ChairSubstituentSemantics, ...]

    def __post_init__(self) -> None:
        if self.conformation not in ("chair_a", "chair_b"):
            raise ValueError("invalid chair conformation")
        if len(self.ring_atom_indices) != 6:
            raise ValueError("chair semantics require six ring atoms")
        _validate_indices(self.ring_atom_indices, name="chair ring indices")
        if len(set(self.ring_atom_indices)) != 6:
            raise ValueError("chair ring atom indices must be unique")
        positions = [item.ring_position for item in self.substituents]
        if len(positions) != len(set(positions)):
            raise ValueError("chair semantic substituent positions must be unique")
        for item in self.substituents:
            if item.ring_atom_index != self.ring_atom_indices[item.ring_position]:
                raise ValueError("chair substituent ring index does not match its position")

    @property
    def referenced_atom_indices(self) -> tuple[int, ...]:
        return self.ring_atom_indices + tuple(
            index
            for substituent in self.substituents
            for index in (substituent.ring_atom_index, substituent.group_atom_index)
        )


ProjectionSemantics: TypeAlias = (
    HaworthProjectionSemantics | FischerProjectionSemantics | ChairProjectionSemantics
)
