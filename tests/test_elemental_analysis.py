"""Tests for formula-string elemental analysis calculations."""

from __future__ import annotations


import pytest


from chemuson.core.elemental_analysis import (
    FormulaError,
    compare_found_calculated,
    elemental_percentages,
    find_solvate_candidates,
    format_formula,
    format_publication_line,
    molecular_weight,
    parse_formula,
)


def assert_composition(actual: dict[str, float], expected: dict[str, float]) -> None:
    assert set(actual) == set(expected)
    for element, count in expected.items():
        assert actual[element] == pytest.approx(count)


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("C6H6", {"C": 6, "H": 6}),
        ("C20H18N4O2", {"C": 20, "H": 18, "N": 4, "O": 2}),
        ("C0.5H1.5", {"C": 0.5, "H": 1.5}),
    ],
)
def test_basic_formula_parsing(formula: str, expected: dict[str, float]) -> None:
    assert_composition(parse_formula(formula), expected)


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("CH2Cl2", {"C": 1, "H": 2, "Cl": 2}),
        ("C2H5Br", {"C": 2, "H": 5, "Br": 1}),
        ("NaCl", {"Na": 1, "Cl": 1}),
        ("FeCl3", {"Fe": 1, "Cl": 3}),
    ],
)
def test_two_letter_elements(formula: str, expected: dict[str, float]) -> None:
    assert_composition(parse_formula(formula), expected)


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("Fe2(SO4)3", {"Fe": 2, "S": 3, "O": 12}),
        ("Cu(NO3)2", {"Cu": 1, "N": 2, "O": 6}),
    ],
)
def test_parentheses(formula: str, expected: dict[str, float]) -> None:
    assert_composition(parse_formula(formula), expected)


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("CuSO4\u00b75H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
        ("CuSO4.H2O", {"Cu": 1, "S": 1, "O": 5, "H": 2}),
        ("C20H18N4O2\u00b70.5H2O", {"C": 20, "H": 19, "N": 4, "O": 2.5}),
        (
            "C20H18N4O2 + 0.25 CH3OH",
            {"C": 20.25, "H": 19, "N": 4, "O": 2.25},
        ),
    ],
)
def test_dot_and_adduct_formulas(formula: str, expected: dict[str, float]) -> None:
    assert_composition(parse_formula(formula), expected)


def test_percent_calculations() -> None:
    percentages = elemental_percentages(parse_formula("C6H6"))

    assert molecular_weight(parse_formula("C6H6")) == pytest.approx(78.11184)
    assert percentages["C"] == pytest.approx(92.26, abs=0.01)
    assert percentages["H"] == pytest.approx(7.74, abs=0.01)
    assert sum(percentages.values()) == pytest.approx(100.0)


def test_found_vs_calculated_comparison() -> None:
    calculated = {"C": 70.12, "H": 5.42, "N": 8.17}
    found = {"C": 69.90, "H": 5.93, "N": 8.14}
    rows = compare_found_calculated(calculated, found, {"C": 0.4, "H": 0.4, "N": 0.4})

    assert [row.element for row in rows] == ["C", "H", "N"]
    assert rows[0].delta == pytest.approx(-0.22)
    assert rows[0].passed is True
    assert rows[1].delta == pytest.approx(0.51)
    assert rows[1].passed is False


def test_publication_line() -> None:
    line = format_publication_line(
        "C20H18N4O2\u00b70.5H2O",
        {"C": 67.13, "H": 5.35, "N": 15.66},
        {"C": 66.98, "H": 5.42, "N": 15.50},
        ["C", "H", "N"],
    )

    assert line == (
        "Anal. Calcd for C20H18N4O2\u00b70.5H2O: "
        "C, 67.13; H, 5.35; N, 15.66. Found: C, 66.98; H, 5.42; N, 15.50."
    )


def test_solvate_finder_ranking_prefers_matching_hydrate() -> None:
    found = elemental_percentages(parse_formula("C20H18N4O2\u00b70.5H2O"))
    candidates = find_solvate_candidates(
        "C20H18N4O2",
        {element: found[element] for element in ("C", "H", "N")},
        ["H2O", "MeOH"],
        equivalent_range=(0.25, 1.0),
        step_size=0.25,
        selected_elements=["C", "H", "N"],
        tolerances=0.4,
        max_components=1,
    )

    assert candidates
    best = candidates[0]
    assert best.formula == "C20H18N4O2\u00b70.5H2O"
    assert best.score == pytest.approx(0.0, abs=1e-9)
    assert best.passed is True


def test_solvate_finder_warns_about_overfitting() -> None:
    found = {"C": 50.0, "H": 5.0}
    candidates = find_solvate_candidates(
        "C10H10",
        found,
        ["H2O", "MeOH"],
        equivalent_range=(0.5, 0.5),
        step_size=0.5,
        selected_elements=["C", "H"],
        max_components=2,
    )

    warned = [candidate for candidate in candidates if len(candidate.solvent_equivalents) == 2]
    assert warned
    assert "Multiple solvent components" in warned[0].warning


@pytest.mark.parametrize(
    ("formula", "message"),
    [
        ("", "empty"),
        ("C6H6)", "closing parenthesis"),
        ("Fe2(SO4", "opening parenthesis"),
        ("Xx2", "Unknown element"),
        ("C0H4", "must be positive"),
        ("C6H6 + ", "Missing formula component"),
    ],
)
def test_invalid_formulas_have_user_friendly_messages(formula: str, message: str) -> None:
    with pytest.raises(FormulaError, match=message):
        parse_formula(formula)


def test_format_formula_normalizes_hill_order() -> None:
    assert format_formula(parse_formula("H2O.C6H6")) == "C6H8O"
