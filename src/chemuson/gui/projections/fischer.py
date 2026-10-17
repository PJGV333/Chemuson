"""Parameterized geometry for Fischer projections."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal

from .common import ProjectionGeometry, projection_atom, projection_bond
from .semantics import FischerCenterSemantics, FischerProjectionSemantics

Side = Literal["left", "right"]


@dataclass(frozen=True, slots=True)
class FischerParameters:
    """Vertical backbone, horizontal crosses, and explicit substituent ordering."""

    vertical_backbone_length: float = 64.0
    horizontal_cross_length: float = 64.0
    canonical_bond_angle_degrees: float = 90.0
    substituent_sides: tuple[Side, ...] = ("right", "left", "right", "right")
    center_configurations: tuple[Literal["R", "S"] | None, ...] = ()
    top_terminal_group: str = "CHO"
    bottom_terminal_group: str = "CH2OH"

    def __post_init__(self) -> None:
        for name, value in (
            ("vertical_backbone_length", self.vertical_backbone_length),
            ("horizontal_cross_length", self.horizontal_cross_length),
            ("canonical_bond_angle_degrees", self.canonical_bond_angle_degrees),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.vertical_backbone_length <= 0.0 or self.horizontal_cross_length <= 0.0:
            raise ValueError("Fischer bond lengths must be positive")
        if not 60.0 <= self.canonical_bond_angle_degrees <= 120.0:
            raise ValueError("canonical_bond_angle_degrees must be between 60 and 120")
        if any(side not in ("left", "right") for side in self.substituent_sides):
            raise ValueError("Fischer substituent sides must be 'left' or 'right'")
        if self.center_configurations and len(self.center_configurations) != len(self.substituent_sides):
            raise ValueError("Fischer center_configurations must match the stereocenter count")
        if any(value not in (None, "R", "S") for value in self.center_configurations):
            raise ValueError("Fischer center configuration must be R, S or undeclared")
        if not self.top_terminal_group.strip() or not self.bottom_terminal_group.strip():
            raise ValueError("Fischer terminal group labels must not be empty")


def build_fischer_geometry(parameters: FischerParameters | None = None) -> ProjectionGeometry:
    """Build a deterministic Fischer carbon chain without stereo metadata."""
    config = parameters or FischerParameters()
    count = len(config.substituent_sides)
    if count < 1:
        raise ValueError("Fischer projection requires at least one stereocenter")
    segment = config.vertical_backbone_length
    total_steps = count + 1
    top_y = -0.5 * total_steps * segment
    atoms = [projection_atom(config.top_terminal_group, 0.0, top_y, explicit=True)]
    centers = []
    for index in range(count):
        y = top_y + (index + 1) * segment
        centers.append(len(atoms))
        atoms.append(projection_atom("C", 0.0, y, explicit=True))
    bottom_index = len(atoms)
    bottom_y = top_y + total_steps * segment
    atoms.append(projection_atom(config.bottom_terminal_group, 0.0, bottom_y, explicit=True))

    bonds = [projection_bond(index, index + 1) for index in range(count + 1)]
    theta = math.radians(90.0 - config.canonical_bond_angle_degrees)
    dx_unit = math.cos(theta)
    dy_unit = math.sin(theta)
    semantic_centers = []
    configurations = config.center_configurations or (None,) * count
    for center_number, (center_index, side, absolute_configuration) in enumerate(
        zip(centers, config.substituent_sides, configurations, strict=True)
    ):
        center_y = atoms[center_index].y
        oh_sign = -1.0 if side == "left" else 1.0
        for label, sign in (("OH", oh_sign), ("H", -oh_sign)):
            substituent_index = len(atoms)
            atoms.append(
                projection_atom(
                    label,
                    sign * dx_unit * config.horizontal_cross_length,
                    center_y + sign * dy_unit * config.horizontal_cross_length,
                    explicit=True,
                )
            )
            bonds.append(projection_bond(center_index, substituent_index))
        first_group_index = count + 2 + center_number * 2
        second_group_index = first_group_index + 1
        left_index, right_index = (
            (first_group_index, second_group_index)
            if side == "left"
            else (second_group_index, first_group_index)
        )
        semantic_centers.append(
            FischerCenterSemantics(
                center_atom_index=center_index,
                left_atom_index=left_index,
                right_atom_index=right_index,
                upper_atom_index=center_index - 1,
                lower_atom_index=center_index + 1,
                absolute_configuration=absolute_configuration,
            )
        )

    return ProjectionGeometry(
        family="fischer",
        bond_length=config.vertical_backbone_length,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        semantics=FischerProjectionSemantics(tuple(semantic_centers)),
    )
