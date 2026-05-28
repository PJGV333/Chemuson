"""Pruebas del predictor espectral base."""

from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.spectroscopy import SpectralPrediction, predict_spectra, register_predictor


def test_predict_spectra_returns_nmr_and_mass_for_ethanol() -> None:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 40.0, 0.0)
    o = graph.add_atom("O", 80.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o.id, order=1)

    prediction = predict_spectra(graph)

    assert prediction.source == "heuristic-v1"
    assert prediction.confidence > 0.0
    assert len(prediction.carbon_nmr) == 2
    assert any(peak.environment == "C unido a heteroátomo" for peak in prediction.carbon_nmr)
    assert any(peak.environment == "alquilo unido a heteroátomo" for peak in prediction.proton_nmr)
    assert all(0.0 < peak.confidence <= 1.0 for peak in prediction.proton_nmr)
    assert prediction.mass_spectrum[0].label == "M+"
    assert prediction.mass_spectrum[0].mz > 45.0


def test_predict_spectra_expands_superatoms_for_calculation() -> None:
    graph = MolGraph()
    root = graph.add_atom("C", 0.0, 0.0)
    phenyl = graph.add_atom("Ph", 40.0, 0.0, is_explicit=True)
    graph.add_bond(root.id, phenyl.id, order=1)

    prediction = predict_spectra(graph)

    assert len(prediction.carbon_nmr) >= 7
    assert any(peak.environment == "aromático" for peak in prediction.carbon_nmr)


def test_predict_spectra_uses_registered_plugin() -> None:
    class DummyPredictor:
        name = "dummy"

        def predict(self, graph: MolGraph) -> SpectralPrediction:
            return SpectralPrediction(source="dummy", message=str(len(graph.atoms)))

    register_predictor(DummyPredictor())
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)

    prediction = predict_spectra(graph, predictor="dummy")

    missing = predict_spectra(graph, predictor="missing")

    assert prediction.source == "dummy"
    assert prediction.message == "1"
    assert missing.message == "predictor_not_found:missing"


def test_predict_spectra_distinguishes_aliphatic_carbon_environments() -> None:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 40.0, 0.0)
    c3 = graph.add_atom("C", 80.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, c3.id, order=1)

    prediction = predict_spectra(graph)
    environments = {peak.environment for peak in prediction.carbon_nmr}

    assert "metilo alifático" in environments
    assert "metileno alifático" in environments


def test_predict_spectra_labels_carboxylic_acid_oh() -> None:
    graph = MolGraph()
    c = graph.add_atom("C", 0.0, 0.0)
    o_double = graph.add_atom("O", 40.0, -40.0, is_explicit=True)
    o_h = graph.add_atom("O", 40.0, 40.0, is_explicit=True, explicit_h=1)
    graph.add_bond(c.id, o_double.id, order=2)
    graph.add_bond(c.id, o_h.id, order=1)

    prediction = predict_spectra(graph)

    assert any(peak.environment == "ácido carboxílico O-H" for peak in prediction.proton_nmr)
