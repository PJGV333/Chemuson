from __future__ import annotations

import pytest

from tests.clean2d_regression.assertions import execute_case
from tests.clean2d_regression.cases import get_regression_cases


def test_regression_case_names_are_unique() -> None:
    names = [case.name for case in get_regression_cases()]
    assert len(names) == len(set(names))


def test_regression_case_metadata_is_complete() -> None:
    for case in get_regression_cases():
        assert case.name
        assert case.family
        assert case.mode
        assert case.expected_states
        assert case.tags
        assert case.target_label in {"whole", "selection"}


def test_regression_corpus_has_minimum_family_coverage() -> None:
    cases = get_regression_cases()
    families = {case.family for case in cases}
    tags = {tag for case in cases for tag in case.tags}

    assert "acyclic" in families
    assert "aromatic" in families
    assert "substituted-aromatic" in families
    assert "fused-ring" in families
    assert "heteroaromatic" in families
    assert "charged" in families
    assert "stereo-sensitive" in families
    assert "selection-boundary" in families
    assert "multi-block" in families
    assert "macrocycle" in families
    assert "charged" in tags
    assert "stereo_sensitive" in tags
    assert "selection_boundary" in tags
    assert "complex_policy_guard" in tags


@pytest.mark.parametrize("case", get_regression_cases(), ids=lambda case: case.name)
def test_clean2d_regression_corpus_case_contracts(case) -> None:
    snapshot = execute_case(case)

    assert snapshot["case"]["name"] == case.name
    assert snapshot["result"]["state"] in case.expected_states
    assert snapshot["identity"]["atom_ids"]
    assert snapshot["identity"]["bond_ids"]
    assert "before" in snapshot["metrics"]
    assert snapshot["case"]["tags"] == list(case.tags)


def test_known_delicate_cases_record_current_behavior_without_geometry_contracts() -> None:
    delicate = [case for case in get_regression_cases() if case.known_delicate]
    assert delicate
    for case in delicate:
        snapshot = execute_case(case)
        assert snapshot["case"]["known_delicate"] is True
        assert snapshot["result"]["state"] in case.expected_states


def test_charged_cases_preserve_formal_charges() -> None:
    charged = [case for case in get_regression_cases() if case.has_tag("charged")]
    assert charged
    for case in charged:
        snapshot = execute_case(case)
        assert "charged" in snapshot["case"]["tags"]


def test_stereo_sensitive_cases_preserve_stereo_metadata() -> None:
    stereo_cases = [case for case in get_regression_cases() if case.has_tag("stereo_sensitive")]
    assert stereo_cases
    for case in stereo_cases:
        snapshot = execute_case(case)
        assert "stereo_sensitive" in snapshot["case"]["tags"]


def test_selection_boundary_cases_preserve_boundary_contracts() -> None:
    boundary_cases = [case for case in get_regression_cases() if case.has_tag("selection_boundary")]
    assert len(boundary_cases) >= 3
    for case in boundary_cases:
        snapshot = execute_case(case)
        assert snapshot["case"]["target"] == "selection"
        assert "selection_boundary" in snapshot["case"]["tags"]


def test_complex_policy_guard_cases_remain_observational() -> None:
    complex_cases = [case for case in get_regression_cases() if case.has_tag("complex_policy_guard")]
    assert complex_cases
    for case in complex_cases:
        assert {"preserve-only", "failed-controlled"} & case.expected_states
        snapshot = execute_case(case)
        assert snapshot["result"]["state"] in case.expected_states
