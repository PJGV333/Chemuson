"""Pruebas de superátomos y abreviaturas colapsadas."""

from __future__ import annotations


from PyQt6.QtWidgets import QApplication


from chemuson.chemio.rdkit_io import (
    SUPPORTED_ABBREVIATION_LABELS,
    abbreviation_expansion_elements,
    expand_abbreviations_for_calculation,
    is_supported_abbreviation_label,
)
from chemuson.core.model import MolGraph
from chemuson.gui.canvas import ChemusonCanvas


def test_supported_abbreviation_metadata_is_case_insensitive() -> None:
    assert "Me" in SUPPORTED_ABBREVIATION_LABELS
    assert "Boc" in SUPPORTED_ABBREVIATION_LABELS
    assert "TBS" in SUPPORTED_ABBREVIATION_LABELS
    assert is_supported_abbreviation_label("tbs")
    assert abbreviation_expansion_elements("TBS") == ("Si", "C", "C", "C", "C", "C", "C")
    assert not is_supported_abbreviation_label("UnsupportedGroup")
    assert abbreviation_expansion_elements("UnsupportedGroup") == ()


def test_tbs_expands_for_calculation_without_mutating_visible_label() -> None:
    graph = MolGraph()
    oxygen = graph.add_atom("O", 0.0, 0.0, is_explicit=True)
    tbs = graph.add_atom("TBS", 40.0, 0.0, is_explicit=True)
    graph.add_bond(oxygen.id, tbs.id, order=1)

    expanded = expand_abbreviations_for_calculation(graph)
    elements = [atom.element for atom in expanded.atoms.values()]

    assert graph.get_atom(tbs.id).element == "TBS"
    assert "TBS" not in elements
    assert elements.count("Si") == 1
    assert elements.count("C") == 6
    assert len(expanded.bonds) == 7


def test_unsupported_abbreviation_degrades_as_collapsed_label() -> None:
    graph = MolGraph()
    label = graph.add_atom("Foo", 0.0, 0.0, is_explicit=True)

    expanded = expand_abbreviations_for_calculation(graph)

    assert expanded is graph
    assert expanded.get_atom(label.id).element == "Foo"


def test_analysis_formula_expands_supported_superatoms() -> None:
    QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    root = canvas.model.add_atom("C", 0.0, 0.0)
    phenyl = canvas.model.add_atom("Ph", 40.0, 0.0, is_explicit=True)
    canvas.model.add_bond(root.id, phenyl.id, order=1)

    text = canvas._analysis_build_text(canvas.model, "formula")

    assert "Chemical Formula: C7" in text
    assert "Ph" not in text
    assert "N/D" not in text


def test_expanded_validation_reports_overbonded_anchor_atom() -> None:
    graph = MolGraph()
    ome = graph.add_atom("OMe", 0.0, 0.0, is_explicit=True)
    for idx in range(3):
        carbon = graph.add_atom("C", 40.0 * (idx + 1), 0.0)
        graph.add_bond(ome.id, carbon.id, order=1)

    expanded = expand_abbreviations_for_calculation(graph)
    issues = expanded.validate_detailed()

    assert ome.id in issues
    assert issues[ome.id].code == "VALENCE_EXCEEDED"


def test_canvas_validation_reports_supported_superatom_anchor() -> None:
    QApplication.instance() or QApplication([])
    canvas = ChemusonCanvas()
    ome = canvas.model.add_atom("OMe", 0.0, 0.0, is_explicit=True)
    for idx in range(3):
        carbon = canvas.model.add_atom("C", 40.0 * (idx + 1), 0.0)
        canvas.model.add_bond(ome.id, carbon.id, order=1)
    canvas._rebuild_items_from_model()

    errors = canvas.validate_structure()

    assert errors == [ome.id]
    assert canvas.current_validation_issues()[ome.id].code == "VALENCE_EXCEEDED"
    assert canvas.atom_items[ome.id]._valence_error is True
