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
        assert case.target_label in {"whole", "selection"}


@pytest.mark.parametrize("case", get_regression_cases(), ids=lambda case: case.name)
def test_clean2d_regression_corpus_case_contracts(case) -> None:
    snapshot = execute_case(case)

    assert snapshot["case"]["name"] == case.name
    assert snapshot["result"]["state"] in case.expected_states
    assert snapshot["identity"]["atom_ids"]
    assert snapshot["identity"]["bond_ids"]
    assert "before" in snapshot["metrics"]


def test_known_delicate_cases_record_current_behavior_without_geometry_contracts() -> None:
    delicate = [case for case in get_regression_cases() if case.known_delicate]
    assert delicate
    for case in delicate:
        snapshot = execute_case(case)
        assert snapshot["case"]["known_delicate"] is True
        assert snapshot["result"]["state"] in case.expected_states
