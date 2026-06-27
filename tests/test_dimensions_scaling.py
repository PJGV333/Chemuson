"""Pruebas para redimensionado global y de selección."""

from __future__ import annotations

import math
from dataclasses import replace

import pytest
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QApplication


from chemuson.chemio.persistence import PersistenceManager
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.items import TextAnnotationItem


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


def _distance(p1: QPointF, p2: QPointF) -> float:
    return math.hypot(p2.x() - p1.x(), p2.y() - p1.y())


def test_scale_current_selection_scales_geometry_and_style_with_undo():
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("N", 100.0, 100.0)
    atom_b = canvas.model.add_atom("N", 140.0, 100.0)
    bond = canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    text_item = TextAnnotationItem("Nota", 200.0, 100.0)
    text_font = text_item.font()
    text_font.setPointSizeF(12.0)
    text_item.setFont(text_font)
    canvas.add_text_item(text_item)

    arrow = canvas.add_arrow_item(QPointF(220.0, 100.0), QPointF(280.0, 100.0), "forward")
    bracket = canvas.add_bracket_item(QRectF(300.0, 80.0, 40.0, 60.0), "[]")

    canvas.atom_items[atom_a.id].setSelected(True)
    canvas.atom_items[atom_b.id].setSelected(True)
    text_item.setSelected(True)
    arrow.setSelected(True)
    bracket.setSelected(True)
    canvas._sync_selection_from_scene()

    old_bond_length = math.hypot(atom_b.x - atom_a.x, atom_b.y - atom_a.y)
    old_text_size = text_item.font().pointSizeF()
    old_arrow_length = _distance(arrow.start_point(), arrow.end_point())
    old_bracket_width = bracket.base_rect().width()
    old_padding = float(bracket._padding)
    default_stroke = float(canvas.drawing_style.stroke_px)

    assert canvas.scale_current_selection(1.5, include_style=True)

    scaled_bond = canvas.model.get_bond(bond.id)
    scaled_atom_a = canvas.model.get_atom(atom_a.id)
    scaled_atom_b = canvas.model.get_atom(atom_b.id)

    assert math.hypot(scaled_atom_b.x - scaled_atom_a.x, scaled_atom_b.y - scaled_atom_a.y) == pytest.approx(
        old_bond_length * 1.5
    )
    assert scaled_atom_a.label_scale == pytest.approx(1.5)
    assert scaled_atom_b.label_scale == pytest.approx(1.5)
    assert scaled_bond.stroke_px == pytest.approx(default_stroke * 1.5)
    assert text_item.font().pointSizeF() == pytest.approx(old_text_size * 1.5)
    assert _distance(arrow.start_point(), arrow.end_point()) == pytest.approx(old_arrow_length * 1.5)
    assert arrow.stroke_px() == pytest.approx(default_stroke * 1.5)
    assert bracket.base_rect().width() == pytest.approx(old_bracket_width * 1.5)
    assert float(bracket._padding) == pytest.approx(old_padding * 1.5)
    assert bracket.stroke_px() == pytest.approx(default_stroke * 1.5)

    canvas.undo_stack.undo()

    restored_bond = canvas.model.get_bond(bond.id)
    restored_atom_a = canvas.model.get_atom(atom_a.id)
    restored_atom_b = canvas.model.get_atom(atom_b.id)
    assert math.hypot(
        restored_atom_b.x - restored_atom_a.x,
        restored_atom_b.y - restored_atom_a.y,
    ) == pytest.approx(old_bond_length)
    assert restored_atom_a.label_scale is None
    assert restored_atom_b.label_scale is None
    assert restored_bond.stroke_px is None
    assert text_item.font().pointSizeF() == pytest.approx(old_text_size)
    assert _distance(arrow.start_point(), arrow.end_point()) == pytest.approx(old_arrow_length)
    assert arrow.stroke_px() is None
    assert bracket.base_rect().width() == pytest.approx(old_bracket_width)
    assert float(bracket._padding) == pytest.approx(old_padding)
    assert bracket.stroke_px() is None


def test_apply_document_dimensions_rescales_existing_content():
    canvas = ChemusonCanvas()
    atom_a = canvas.model.add_atom("N", 120.0, 120.0)
    atom_b = canvas.model.add_atom("N", 160.0, 120.0)
    canvas.model.add_bond(atom_a.id, atom_b.id)
    canvas._rebuild_items_from_model()

    text_item = TextAnnotationItem("Escala", 210.0, 120.0)
    text_font = text_item.font()
    text_font.setPointSizeF(10.0)
    text_item.setFont(text_font)
    canvas.add_text_item(text_item)
    arrow = canvas.add_arrow_item(QPointF(240.0, 120.0), QPointF(300.0, 120.0), "forward")

    old_bond_length = math.hypot(atom_b.x - atom_a.x, atom_b.y - atom_a.y)
    old_text_size = text_item.font().pointSizeF()
    old_arrow_length = _distance(arrow.start_point(), arrow.end_point())

    factor = 1.5
    style = replace(
        canvas.drawing_style,
        bond_length_px=canvas.drawing_style.bond_length_px * factor,
        stroke_px=canvas.drawing_style.stroke_px * factor,
        double_offset_px=canvas.drawing_style.double_offset_px * factor,
        wedge_width_px=canvas.drawing_style.wedge_width_px * factor,
        hash_stroke_px=canvas.drawing_style.hash_stroke_px * factor,
    )

    canvas.apply_document_dimensions(
        style=style,
        label_font_size=canvas.state.label_font_size * factor,
        numbering_font_size=canvas.state.numbering_font_size * factor,
        scale_existing=True,
        scale_factor=factor,
    )

    scaled_atom_a = canvas.model.get_atom(atom_a.id)
    scaled_atom_b = canvas.model.get_atom(atom_b.id)
    assert math.hypot(scaled_atom_b.x - scaled_atom_a.x, scaled_atom_b.y - scaled_atom_a.y) == pytest.approx(
        old_bond_length * factor
    )
    assert canvas.state.bond_length == pytest.approx(style.bond_length_px)
    assert canvas.state.label_font_size == pytest.approx(11.0 * factor)
    assert text_item.font().pointSizeF() == pytest.approx(old_text_size * factor)
    assert _distance(arrow.start_point(), arrow.end_point()) == pytest.approx(old_arrow_length * factor)
    assert arrow.stroke_px() is None


def test_persistence_roundtrip_keeps_local_label_scale_bracket_stroke_and_text_width():
    canvas = ChemusonCanvas()
    atom = canvas.model.add_atom("N", 100.0, 100.0, label_scale=1.25)
    canvas._rebuild_items_from_model()

    text_item = TextAnnotationItem("Bloque", 40.0, 50.0)
    text_item.setTextWidth(120.0)
    canvas.add_text_item(text_item)
    bracket = canvas.add_bracket_item(QRectF(10.0, 10.0, 40.0, 50.0), "[]", stroke_px=4.2)

    data = PersistenceManager.save_to_dict(canvas)

    restored = ChemusonCanvas()
    PersistenceManager.load_from_dict(data, restored)

    restored_atom = restored.model.get_atom(atom.id)
    assert restored_atom.label_scale == pytest.approx(1.25)
    assert len(restored.bracket_items) == 1
    assert restored.bracket_items[0].stroke_px() == pytest.approx(4.2)
    restored_texts = restored._all_scalable_text_items()
    assert len(restored_texts) == 1
    assert restored_texts[0].textWidth() == pytest.approx(120.0)
