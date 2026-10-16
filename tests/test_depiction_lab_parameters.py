from __future__ import annotations

import json

import pytest

from chemuson.gui.style import CHEMDOODLE_LIKE
from tools.visual.depiction_lab_app.parameters import (
    PARAMETER_FIELDS,
    DepictionStyleParameters,
    ParameterValidationError,
    REFERENCE_BOND_LENGTH,
)
from tools.visual.depiction_lab_app.profiles import (
    ProfileError,
    compare_profiles,
    decode_profile_text,
    encode_profile,
    load_profile,
    profile_directory,
    profile_payload,
    save_profile,
)


def test_current_profile_maps_exactly_to_existing_drawing_style() -> None:
    parameters = DepictionStyleParameters.current()
    assert set(parameters.to_mapping()) == set(PARAMETER_FIELDS)
    assert parameters.to_drawing_style() == CHEMDOODLE_LIKE
    assert REFERENCE_BOND_LENGTH == CHEMDOODLE_LIKE.bond_length_px
    assert all(
        getattr(parameters, name) is None
        for name in (
            "triple_bond_offset_ratio",
            "bold_double_bond_offset_ratio",
            "solid_wedge_tip_width_ratio",
            "hashed_wedge_terminal_shift_ratio",
            "label_clearance_ratio",
        )
    )


def test_profile_json_roundtrip_is_stable_and_profile_diff_is_explicit() -> None:
    current = DepictionStyleParameters.current()
    values = current.to_mapping()
    values["solid_wedge_width_ratio"] = 0.31
    values["hashed_wedge_count"] = 8
    values["label_clearance_ratio"] = 0.035
    working = DepictionStyleParameters.from_mapping(values)
    encoded = encode_profile(working, name="Publication candidate")
    name, loaded = decode_profile_text(encoded)

    assert name == "Publication candidate"
    assert loaded == working
    assert encode_profile(loaded, name=name) == encoded
    assert set(compare_profiles(current, working)) == {
        "solid_wedge_width_ratio",
        "hashed_wedge_count",
        "label_clearance_ratio",
    }


def test_profile_rejects_invalid_json_schema_types_ranges_and_unknown_fields() -> None:
    valid = profile_payload(DepictionStyleParameters.current())
    bad_cases = []
    bad = dict(valid)
    bad["unknown"] = 1
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["version"] = 90
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["not_a_renderer_control"] = 1.0
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["bond_stroke_ratio"] = True
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["hashed_wedge_count"] = 6.0
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["solid_wedge_width_ratio"] = 9.0
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["hashed_wedge_min_width_ratio"] = 0.30
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["parameters"]["cap_style"] = []
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["name"] = "  "
    bad_cases.append(bad)
    bad = json.loads(json.dumps(valid))
    bad["based_on"] = ""
    bad_cases.append(bad)

    for payload in bad_cases:
        with pytest.raises(ProfileError):
            decode_profile_text(json.dumps(payload))
    with pytest.raises(ProfileError, match="non-finite"):
        decode_profile_text('{"schema":"chemuson.depiction-style","version":1,"parameters":{"x":NaN}}')
    with pytest.raises(ProfileError, match="invalid profile JSON"):
        decode_profile_text("{")
    with pytest.raises(ProfileError, match="duplicate JSON field"):
        decode_profile_text(
            '{"schema":"chemuson.depiction-style","schema":"chemuson.depiction-style"}'
        )


def test_parameter_model_rejects_missing_and_cross_parameter_invariants() -> None:
    values = DepictionStyleParameters.current().to_mapping()
    with pytest.raises(ParameterValidationError, match="missing"):
        DepictionStyleParameters.from_mapping({"bond_stroke_ratio": 0.065})
    values["hashed_wedge_min_width_ratio"] = 0.14
    values["hashed_wedge_max_width_ratio"] = 0.10
    with pytest.raises(ParameterValidationError, match="minimum width"):
        DepictionStyleParameters.from_mapping(values)


def test_version_zero_profile_uses_explicit_compatibility_migration() -> None:
    old = {
        "schema": "chemuson.depiction-style",
        "version": 0,
        "name": "Old candidate",
        "based_on": "current-renderer",
        "parameters": {"wedge_width_ratio": 0.30, "hash_count": 7},
    }
    name, migrated = decode_profile_text(json.dumps(old))
    assert name == "Old candidate"
    assert migrated.solid_wedge_width_ratio == 0.30
    assert migrated.hashed_wedge_count == 7
    assert migrated.double_bond_offset_ratio == DepictionStyleParameters.current().double_bond_offset_ratio
    old["unexpected"] = 4
    with pytest.raises(ProfileError, match="unknown version 0"):
        decode_profile_text(json.dumps(old))


def test_profile_save_reload_and_lab_owned_directory(tmp_path) -> None:
    values = DepictionStyleParameters.current().to_mapping()
    values["double_bond_offset_ratio"] = 0.16
    profile = DepictionStyleParameters.from_mapping(values)
    path = save_profile(tmp_path / "profiles" / "candidate.json", profile, name="Candidate")
    name, loaded = load_profile(path)

    assert name == "Candidate"
    assert loaded == profile
    assert profile_directory().parts[-2:] == ("depiction-lab", "profiles")
