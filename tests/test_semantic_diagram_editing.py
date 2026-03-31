"""Pruebas de edición/reconstrucción para diagramas semánticos."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import EditSemanticDiagramCommand
from chemuson.gui.energy_diagrams import (
    build_atomic_species_diagram,
    build_atomic_subshell_diagram,
)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def test_edit_semantic_diagram_command_preserves_transform_and_updates_summary() -> None:
    canvas = ChemusonCanvas()
    try:
        inserted = canvas.insert_semantic_diagram(
            build_atomic_subshell_diagram(8),
            QPointF(220.0, 180.0),
        )
        item = inserted[0]
        item.set_display_rect(QRectF(80.0, 95.0, 220.0, 160.0))
        item.setRotation(17.5)
        item.setZValue(81.0)
        canvas.set_graphics_item_opacity(item, 0.55)
        before_rect = item.display_rect()
        before_rotation = item.rotation()
        before_z = item.zValue()
        before_opacity = canvas.item_raw_opacity(item)

        before_payload = item.to_json()
        before_payload["opacity"] = before_opacity
        after_payload = dict(before_payload)
        after_payload["semantic_diagram"] = build_atomic_subshell_diagram(7).to_json_dict()

        command = EditSemanticDiagramCommand(canvas, item, before_payload, after_payload)
        canvas.undo_stack.push(command)

        assert item.display_rect() == before_rect
        assert item.rotation() == before_rotation
        assert item.zValue() == before_z
        assert canvas.item_raw_opacity(item) == before_opacity
        assert item.semantic_diagram.metadata["summary_lines"][0] == "1s2 2s2 2p3"

        canvas.undo_stack.undo()

        assert item.display_rect() == before_rect
        assert item.rotation() == before_rotation
        assert item.zValue() == before_z
        assert canvas.item_raw_opacity(item) == before_opacity
        assert item.semantic_diagram.metadata["summary_lines"][0] == "1s2 2s2 2p4"
    finally:
        canvas.close()


def test_semantic_diagram_builder_metadata_roundtrips_in_persistence() -> None:
    canvas = ChemusonCanvas()
    restored = None
    try:
        canvas.insert_semantic_diagram(
            build_atomic_species_diagram("Cr"),
            QPointF(250.0, 200.0),
        )
        data = canvas.get_persistence_data()
        restored = ChemusonCanvas()
        restored.load_persistence_data(data)

        assert len(restored.semantic_diagram_items) == 1
        metadata = restored.semantic_diagram_items[0].semantic_diagram.metadata
        assert metadata["builder"]["name"] == "build_atomic_species_diagram"
        assert metadata["builder"]["params"]["symbol"] == "Cr"
        assert metadata["configuration_string"].endswith("3d5 4s1")
    finally:
        if restored is not None:
            restored.close()
        canvas.close()


def test_semantic_diagram_summary_can_be_hidden_with_undo() -> None:
    canvas = ChemusonCanvas()
    try:
        inserted = canvas.insert_semantic_diagram(
            build_atomic_subshell_diagram(8),
            QPointF(240.0, 180.0),
        )
        item = inserted[0]

        assert item.semantic_diagram.metadata["show_summary"] is True

        changed = canvas._set_semantic_diagram_summary_visible(item, False)

        assert changed is True
        assert item.semantic_diagram.metadata["show_summary"] is False

        canvas.undo_stack.undo()

        assert item.semantic_diagram.metadata["show_summary"] is True
    finally:
        canvas.close()
