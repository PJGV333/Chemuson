"""Opt-in, reusable named profiles for the isolated depiction lab."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

from .parameters import DepictionStyleParameters

EDITORIAL_AB3_NAME = "ChemUSON Editorial AB3"
EDITORIAL_AB3_BASED_ON = (
    "owner-approved valid cases A1–B3; B4 excluded; experimental only; "
    "not a production default"
)


@dataclass(frozen=True, slots=True)
class NamedDepictionPreset:
    name: str
    parameters: DepictionStyleParameters
    based_on: str


EDITORIAL_AB3 = NamedDepictionPreset(
    name=EDITORIAL_AB3_NAME,
    parameters=DepictionStyleParameters.from_mapping(
        {
            "bond_stroke_ratio": 0.075,
            "double_bond_offset_ratio": 0.150,
            "double_bond_inner_trim_ratio": 0.175,
            "triple_bond_offset_ratio": None,
            "bold_double_bond_offset_ratio": None,
            "solid_wedge_width_ratio": 0.450,
            "solid_wedge_tip_width_ratio": None,
            "hashed_wedge_count": 9,
            "hashed_wedge_min_width_ratio": 0.050,
            "hashed_wedge_max_width_ratio": 0.370,
            "hashed_wedge_stroke_ratio": 0.055,
            "hashed_wedge_terminal_shift_ratio": None,
            "label_clearance_ratio": 0.020,
            "cap_style": "round",
            "join_style": "round",
        }
    ),
    based_on=EDITORIAL_AB3_BASED_ON,
)

NAMED_PRESETS: Mapping[str, NamedDepictionPreset] = MappingProxyType(
    {EDITORIAL_AB3_NAME: EDITORIAL_AB3}
)


def get_named_preset(name: str) -> NamedDepictionPreset:
    """Return a known lab preset or fail with available names."""
    try:
        return NAMED_PRESETS[name]
    except KeyError as exc:
        available = ", ".join(NAMED_PRESETS) or "none"
        raise KeyError(f"unknown depiction preset {name!r}; available: {available}") from exc
