"""Regresiones de la vista GUI de espectros."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.docks import ChemicalPropertiesDock, SpectroscopyDock
from chemuson.gui.main_window import ChemusonWindow
from chemuson.spectroscopy import CarbonNmrPeak, MassPeak, ProtonNmrPeak, SpectralPrediction


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_spectroscopy_dock_populates_peak_tables_and_emits_atom() -> None:
    dock = SpectroscopyDock()
    selected: list[int] = []
    dock.peak_atom_selected.connect(selected.append)

    dock.update_prediction(
        SpectralPrediction(
            proton_nmr=[ProtonNmrPeak(7, 1.23, 3, "alquilo")],
            carbon_nmr=[CarbonNmrPeak(7, 25.0, "alifático")],
            mass_spectrum=[MassPeak(46.0419, 100.0, "M+")],
        )
    )
    dock.proton_table.selectRow(0)

    assert dock.tabs.isVisible() or dock.proton_table.rowCount() == 1
    assert dock.proton_table.item(0, 0).text() == "1.23"
    assert dock.proton_table.item(0, 3).text() == "0.45"
    assert dock.carbon_table.item(0, 1).text() == "alifático"
    assert dock.mass_table.item(0, 2).text() == "M+"
    assert selected == [7]


def test_spectroscopy_dock_copies_active_table_tsv() -> None:
    dock = SpectroscopyDock()
    dock.update_prediction(
        SpectralPrediction(
            proton_nmr=[ProtonNmrPeak(3, 2.50, 1, "alcohol/fenol O-H intercambiable", 0.42)],
            source="heuristic-v1",
            confidence=0.42,
        )
    )

    text = dock.copy_current_table_tsv()

    assert "δ ppm\tInt.\tEntorno\tConf.\tÁtomo" in text
    assert "2.50\t1H\talcohol/fenol O-H intercambiable\t0.42\t3" in text
    assert QApplication.clipboard().text() == text


def test_spectroscopy_dock_exports_active_table_tsv(tmp_path) -> None:
    dock = SpectroscopyDock()
    dock.update_prediction(
        SpectralPrediction(
            carbon_nmr=[CarbonNmrPeak(2, 171.0, "carbonilo", 0.62)],
            source="heuristic-v1",
            confidence=0.62,
        )
    )
    dock.tabs.setCurrentWidget(dock.carbon_table)

    out = tmp_path / "spectrum.tsv"
    path = dock.export_current_table_tsv(str(out))

    assert path == str(out)
    assert "171.0\tcarbonilo\t0.62\t2" in out.read_text(encoding="utf-8")


def test_chemical_properties_dock_copies_tsv() -> None:
    dock = ChemicalPropertiesDock()
    dock.update_properties([("Fórmula", "C2H6O"), ("logP", "-0.14")])

    text = dock.copy_properties_tsv()

    assert text.splitlines()[0] == "Propiedad\tValor"
    assert "Fórmula\tC2H6O" in text
    assert QApplication.clipboard().text() == text


def test_spectrum_peak_selection_selects_canvas_atom() -> None:
    window = ChemusonWindow()
    atom = window.canvas.model.add_atom("C", 100.0, 100.0)
    window.canvas.add_atom_item(atom)

    window._select_atom_from_spectrum(atom.id)

    assert window.canvas.state.selected_atoms == {atom.id}
    assert window.canvas.atom_items[atom.id].isSelected()
