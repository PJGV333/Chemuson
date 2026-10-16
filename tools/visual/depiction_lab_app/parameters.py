"""Typed and validated visual-lab depiction parameters."""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, Mapping

from PyQt6.QtCore import Qt

from chemuson.gui.style import CHEMDOODLE_LIKE, DrawingStyle

SCHEMA = "chemuson.depiction-style"
SCHEMA_VERSION = 1
REFERENCE_BOND_LENGTH = float(CHEMDOODLE_LIKE.bond_length_px)

CAP_STYLES = {
    "flat": Qt.PenCapStyle.FlatCap,
    "square": Qt.PenCapStyle.SquareCap,
    "round": Qt.PenCapStyle.RoundCap,
}
JOIN_STYLES = {
    "miter": Qt.PenJoinStyle.MiterJoin,
    "bevel": Qt.PenJoinStyle.BevelJoin,
    "round": Qt.PenJoinStyle.RoundJoin,
}

# Public profile keys are intentionally explicit and stable.
PARAMETER_FIELDS = (
    "bond_stroke_ratio",
    "double_bond_offset_ratio",
    "double_bond_inner_trim_ratio",
    "triple_bond_offset_ratio",
    "bold_double_bond_offset_ratio",
    "solid_wedge_width_ratio",
    "solid_wedge_tip_width_ratio",
    "hashed_wedge_count",
    "hashed_wedge_min_width_ratio",
    "hashed_wedge_max_width_ratio",
    "hashed_wedge_stroke_ratio",
    "hashed_wedge_terminal_shift_ratio",
    "label_clearance_ratio",
    "cap_style",
    "join_style",
)

RANGES: dict[str, tuple[float, float]] = {
    "bond_stroke_ratio": (0.02, 0.16),
    "double_bond_offset_ratio": (0.03, 0.28),
    "double_bond_inner_trim_ratio": (0.0, 0.50),
    "triple_bond_offset_ratio": (0.03, 0.45),
    "bold_double_bond_offset_ratio": (0.03, 0.45),
    "solid_wedge_width_ratio": (0.10, 0.50),
    "solid_wedge_tip_width_ratio": (0.0, 0.45),
    "hashed_wedge_min_width_ratio": (0.01, 0.15),
    "hashed_wedge_max_width_ratio": (0.05, 0.40),
    "hashed_wedge_stroke_ratio": (0.01, 0.10),
    "hashed_wedge_terminal_shift_ratio": (0.0, 0.50),
    "label_clearance_ratio": (0.005, 0.12),
}
OPTIONAL_FIELDS = frozenset(
    {
        "triple_bond_offset_ratio",
        "bold_double_bond_offset_ratio",
        "solid_wedge_tip_width_ratio",
        "hashed_wedge_terminal_shift_ratio",
        "label_clearance_ratio",
    }
)


class ParameterValidationError(ValueError):
    """A profile contains an unknown, malformed, or unsafe value."""


@dataclass(frozen=True)
class DepictionStyleParameters:
    """Presentation-only values; lengths are ratios of the reference bond."""

    bond_stroke_ratio: float = 0.065
    double_bond_offset_ratio: float = 0.105
    double_bond_inner_trim_ratio: float = 0.175
    triple_bond_offset_ratio: float | None = None
    bold_double_bond_offset_ratio: float | None = None
    solid_wedge_width_ratio: float = 0.225
    solid_wedge_tip_width_ratio: float | None = None
    hashed_wedge_count: int = 6
    hashed_wedge_min_width_ratio: float = 0.05
    hashed_wedge_max_width_ratio: float = 0.25
    hashed_wedge_stroke_ratio: float = 0.055
    hashed_wedge_terminal_shift_ratio: float | None = None
    label_clearance_ratio: float | None = None
    cap_style: str = "round"
    join_style: str = "round"

    @classmethod
    def current(cls) -> "DepictionStyleParameters":
        """Return editable data equal to current renderer values."""
        return cls()

    def to_mapping(self) -> dict[str, Any]:
        return {name: getattr(self, name) for name in PARAMETER_FIELDS}

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> "DepictionStyleParameters":
        if not isinstance(payload, Mapping):
            raise ParameterValidationError("parameters must be a JSON object")
        unknown = sorted(set(payload) - set(PARAMETER_FIELDS))
        missing = sorted(set(PARAMETER_FIELDS) - set(payload))
        if unknown:
            raise ParameterValidationError(f"unknown parameter(s): {', '.join(unknown)}")
        if missing:
            raise ParameterValidationError(f"missing parameter(s): {', '.join(missing)}")

        values = dict(payload)
        for name, (minimum, maximum) in RANGES.items():
            value = values[name]
            if value is None and name in OPTIONAL_FIELDS:
                continue
            if isinstance(value, bool) or not isinstance(value, (float, int)):
                raise ParameterValidationError(f"{name} must be a finite number")
            number = float(value)
            if not math.isfinite(number):
                raise ParameterValidationError(f"{name} must be a finite number")
            if not minimum <= number <= maximum:
                raise ParameterValidationError(
                    f"{name} must be between {minimum:g} and {maximum:g}"
                )
            values[name] = number

        count = values["hashed_wedge_count"]
        if isinstance(count, bool) or not isinstance(count, int):
            raise ParameterValidationError("hashed_wedge_count must be an integer")
        if not 3 <= count <= 12:
            raise ParameterValidationError("hashed_wedge_count must be between 3 and 12")
        if values["hashed_wedge_min_width_ratio"] > values["hashed_wedge_max_width_ratio"]:
            raise ParameterValidationError(
                "hashed wedge minimum width must not exceed maximum width"
            )
        if not isinstance(values["cap_style"], str) or values["cap_style"] not in CAP_STYLES:
            raise ParameterValidationError("cap_style must be flat, square, or round")
        if not isinstance(values["join_style"], str) or values["join_style"] not in JOIN_STYLES:
            raise ParameterValidationError("join_style must be miter, bevel, or round")
        try:
            return cls(**values)
        except TypeError as exc:
            raise ParameterValidationError(str(exc)) from exc

    def to_drawing_style(self, base: DrawingStyle = CHEMDOODLE_LIKE) -> DrawingStyle:
        """Convert profile ratios to the production renderer's style object."""
        length = REFERENCE_BOND_LENGTH
        def optional_px(value: float | None) -> float | None:
            return None if value is None else float(value) * length

        return replace(
            base,
            stroke_px=self.bond_stroke_ratio * length,
            double_offset_px=self.double_bond_offset_ratio * length,
            inner_trim_px=self.double_bond_inner_trim_ratio * length,
            triple_bond_offset_px=optional_px(self.triple_bond_offset_ratio),
            bold_double_offset_px=optional_px(self.bold_double_bond_offset_ratio),
            wedge_width_px=self.solid_wedge_width_ratio * length,
            wedge_tip_width_ratio=self.solid_wedge_tip_width_ratio,
            hash_count=self.hashed_wedge_count,
            hash_min_px=self.hashed_wedge_min_width_ratio * length,
            hash_max_px=self.hashed_wedge_max_width_ratio * length,
            hash_stroke_px=self.hashed_wedge_stroke_ratio * length,
            hashed_terminal_shift_ratio=self.hashed_wedge_terminal_shift_ratio,
            label_clearance_ratio=self.label_clearance_ratio,
            cap_style=CAP_STYLES[self.cap_style],
            join_style=JOIN_STYLES[self.join_style],
        )

    def changed_values(self, other: "DepictionStyleParameters") -> dict[str, dict[str, Any]]:
        """Return a stable profile diff as ``{name: {current, working}}``."""
        left, right = self.to_mapping(), other.to_mapping()
        return {
            name: {"current": left[name], "working": right[name]}
            for name in PARAMETER_FIELDS
            if left[name] != right[name]
        }


def validate_parameters(payload: Mapping[str, Any]) -> DepictionStyleParameters:
    """Named validation entrypoint for UI and profile loaders."""
    return DepictionStyleParameters.from_mapping(payload)
