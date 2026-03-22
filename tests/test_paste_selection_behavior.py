"""Regresiones para selección automática tras pegar."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import TextAnnotationItem


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_paste_selection_payload_selects_newly_pasted_structure_and_text():
    """El pegado interno debe dejar seleccionados los ítems recién creados."""
    source = ChemusonCanvas()
    a1 = source.model.add_atom("N", 100.0, 100.0)
    a2 = source.model.add_atom("O", 140.0, 100.0)
    source.model.add_bond(a1.id, a2.id)
    source._rebuild_items_from_model()

    text_item = TextAnnotationItem("Nota", 170.0, 100.0)
    source.add_text_item(text_item)

    source.atom_items[a1.id].setSelected(True)
    source.atom_items[a2.id].setSelected(True)
    text_item.setSelected(True)
    source._sync_selection_from_scene()

    payload = source._build_selection_payload()
    assert payload is not None

    target = ChemusonCanvas()
    existing = target.model.add_atom("Cl", 40.0, 40.0)
    target._rebuild_items_from_model()
    target.atom_items[existing.id].setSelected(True)
    target._sync_selection_from_scene()

    target._last_scene_pos = QPointF(260.0, 260.0)
    target._paste_selection_payload(payload)

    assert existing.id not in target.state.selected_atoms
    assert len(target.state.selected_atoms) == 2
    selected_items = target.scene.selectedItems()
    assert len(selected_items) == 3
    assert sum(1 for item in selected_items if isinstance(item, TextAnnotationItem)) == 1


def test_insert_molgraph_with_selection_flag_selects_inserted_atoms():
    """El pegado de molfile/SMILES debe seleccionar la estructura insertada."""
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0)
    a2 = graph.add_atom("C", 40.0, 0.0)
    graph.add_bond(a1.id, a2.id)

    canvas = ChemusonCanvas()
    existing = canvas.model.add_atom("Br", 40.0, 40.0)
    canvas._rebuild_items_from_model()
    canvas.atom_items[existing.id].setSelected(True)
    canvas._sync_selection_from_scene()

    canvas._insert_molgraph(graph, select_inserted=True)

    assert existing.id not in canvas.state.selected_atoms
    assert len(canvas.state.selected_atoms) == 2
    assert len(canvas.scene.selectedItems()) == 2
