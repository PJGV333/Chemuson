"""Regresiones de la vista GUI de espectros."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.docks import SpectroscopyDock
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
    assert dock.carbon_table.item(0, 1).text() == "alifático"
    assert dock.mass_table.item(0, 2).text() == "M+"
    assert selected == [7]


def test_spectrum_peak_selection_selects_canvas_atom() -> None:
    window = ChemusonWindow()
    atom = window.canvas.model.add_atom("C", 100.0, 100.0)
    window.canvas.add_atom_item(atom)

    window._select_atom_from_spectrum(atom.id)

    assert window.canvas.state.selected_atoms == {atom.id}
    assert window.canvas.atom_items[atom.id].isSelected()
