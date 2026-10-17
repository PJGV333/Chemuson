from __future__ import annotations

import json

from chemuson.gui.style import CHEMDOODLE_LIKE
from tools.visual.depiction_lab_app.named_presets import (
    EDITORIAL_AB3,
    EDITORIAL_AB3_BASED_ON,
    EDITORIAL_AB3_NAME,
    get_named_preset,
)
from tools.visual.depiction_lab_app.parameters import DepictionStyleParameters
from tools.visual.depiction_lab_app.profiles import decode_profile_text, encode_profile

EXPECTED_EDITORIAL_AB3 = {
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


def test_named_editorial_ab3_contains_exact_owner_approved_values() -> None:
    preset = get_named_preset("ChemUSON Editorial AB3")
    assert preset is EDITORIAL_AB3
    assert preset.name == EDITORIAL_AB3_NAME
    assert preset.based_on == EDITORIAL_AB3_BASED_ON
    assert preset.parameters.to_mapping() == EXPECTED_EDITORIAL_AB3


def test_preset_is_opt_in_and_does_not_mutate_production_or_current_defaults() -> None:
    production_style = CHEMDOODLE_LIKE
    current = DepictionStyleParameters.current()
    current_style = current.to_drawing_style()

    assert current_style == production_style
    assert CHEMDOODLE_LIKE is production_style
    assert current == DepictionStyleParameters.current()
    assert EDITORIAL_AB3.parameters.to_drawing_style() != production_style
    assert EDITORIAL_AB3.parameters.to_mapping() == EXPECTED_EDITORIAL_AB3
    assert current.to_mapping() != EXPECTED_EDITORIAL_AB3


def test_named_preset_profile_roundtrips_with_experimental_provenance() -> None:
    encoded = encode_profile(
        EDITORIAL_AB3.parameters,
        name=EDITORIAL_AB3.name,
        based_on=EDITORIAL_AB3.based_on,
    )
    name, restored = decode_profile_text(encoded)
    payload = json.loads(encoded)

    assert name == EDITORIAL_AB3_NAME
    assert restored == EDITORIAL_AB3.parameters
    assert payload["based_on"] == EDITORIAL_AB3_BASED_ON
    assert "not a production default" in payload["based_on"]
