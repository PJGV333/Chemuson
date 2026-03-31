"""Pruebas para diagramas de niveles de energía persistentes."""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.energy_diagrams import (
    ENERGY_DIAGRAM_MENU_ORDER,
    energy_diagram_display_name,
    energy_diagram_tool_id,
)
from chemuson.gui.main_window import ChemusonWindow


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def test_energy_diagram_tool_id_maps_to_generic_canvas_tool() -> None:
    canvas = ChemusonCanvas()
    try:
        canvas.set_current_tool(energy_diagram_tool_id("hybrid_sp3"))

        assert canvas.current_tool == "tool_energy_diagram"
        assert canvas.state.active_tool == "tool_energy_diagram"
        assert canvas.state.active_energy_diagram_kind == "hybrid_sp3"
    finally:
        canvas.close()


def test_inserted_energy_diagram_is_persistent_and_undoable() -> None:
    canvas = ChemusonCanvas()
    restored = None
    try:
        canvas.state.active_energy_diagram_kind = "sublevel_p"
        center = QPointF(180.0, 140.0)

        item = canvas._insert_energy_diagram_item(center)

        assert item is not None
        assert len(canvas.energy_diagram_items) == 1
        assert item.kind() == "sublevel_p"
        assert item.effective_style()["fill_color"] == "#FFFFFF"
        assert item.effective_style()["box_base_visible"] is False
        assert item.display_rect().center().x() == pytest.approx(center.x())
        assert item.display_rect().center().y() == pytest.approx(center.y())

        item.set_label("2p")
        item.set_label_side("left")
        item.set_occupancies(["up", "pair", "empty"])
        item.set_style_payload(
            {
                "stroke_color": "#112233",
                "fill_color": "#ddeeff",
                "label_color": "#2255ff",
                "arrow_up_color": "#111111",
                "arrow_down_color": "#cc0000",
                "fill_visible": False,
                "box_base_visible": True,
            }
        )
        canvas.set_graphics_item_opacity(item, 0.42)

        data = canvas.get_persistence_data()
        restored = ChemusonCanvas()
        restored.load_persistence_data(data)

        assert len(restored.energy_diagram_items) == 1
        restored_item = restored.energy_diagram_items[0]
        assert restored_item.kind() == "sublevel_p"
        assert restored_item.label() == "2p"
        assert restored_item.label_side() == "left"
        assert restored_item.occupancies() == ("up", "pair", "empty")
        assert restored_item.style_payload() == {
            "stroke_color": "#112233",
            "fill_color": "#ddeeff",
            "label_color": "#2255ff",
            "arrow_up_color": "#111111",
            "arrow_down_color": "#cc0000",
            "fill_visible": False,
            "box_base_visible": True,
        }
        assert restored.effective_item_opacity(restored_item) == pytest.approx(0.42)

        canvas.undo_stack.undo()
        assert len(canvas.energy_diagram_items) == 0

        canvas.undo_stack.redo()
        assert len(canvas.energy_diagram_items) == 1
        assert canvas.energy_diagram_items[0].kind() == "sublevel_p"
    finally:
        if restored is not None:
            restored.close()
        canvas.close()


def test_custom_energy_row_and_legacy_presets_are_supported() -> None:
    canvas = ChemusonCanvas()
    try:
        custom = canvas._insert_energy_diagram_item(QPointF(90.0, 180.0), kind="custom_level")
        aufbau = canvas._insert_energy_diagram_item(QPointF(180.0, 180.0), kind="levels_aufbau")
        mo = canvas._insert_energy_diagram_item(QPointF(420.0, 180.0), kind="mo_2s_2p")

        assert custom is not None
        assert aufbau is not None
        assert mo is not None
        assert custom.kind() == "custom_level"
        assert custom.box_count() == 1
        assert custom.effective_style()["fill_color"] == "#FFFFFF"
        assert custom.effective_style()["box_base_visible"] is False
        assert "custom_level" in ENERGY_DIAGRAM_MENU_ORDER
        assert "levels_aufbau" not in ENERGY_DIAGRAM_MENU_ORDER
        assert "mo_2s_2p" not in ENERGY_DIAGRAM_MENU_ORDER
        assert aufbau.kind() == "levels_aufbau"
        assert aufbau.box_count() == 56
        assert aufbau.effective_style()["fill_color"] == "#FFFFFF"
        assert aufbau.effective_style()["box_stroke_visible"] is True
        assert mo.kind() == "mo_2s_2p"
        assert mo.box_count() == 16
        assert mo.effective_style()["box_stroke_visible"] is False
        assert mo.effective_style()["fill_visible"] is False
    finally:
        canvas.close()


def test_energy_diagram_toolbar_selection_updates_active_canvas() -> None:
    window = ChemusonWindow()
    try:
        tool_id = energy_diagram_tool_id("hybrid_sp2")
        window.symbols_toolbar._select_energy_diagram_tool(tool_id)

        assert window._current_tool_id == tool_id
        assert window.canvas.current_tool == "tool_energy_diagram"
        assert window.canvas.state.active_energy_diagram_kind == "hybrid_sp2"
        assert window.symbols_toolbar.energy_diagram_action.isChecked() is True
        assert (
            window.symbols_toolbar.energy_diagram_action.toolTip()
            == energy_diagram_display_name("hybrid_sp2")
        )
    finally:
        window.canvas.undo_stack.clear()
        window.canvas.undo_stack.setClean()
        window.close()


def test_energy_diagram_copy_payload_roundtrip() -> None:
    source = ChemusonCanvas()
    target = ChemusonCanvas()
    try:
        item = source._insert_energy_diagram_item(QPointF(150.0, 120.0), kind="custom_level")

        assert item is not None
        item.set_slot_count(5)
        item.set_label("4d")
        item.set_label_side("right")
        item.set_occupancies(["pair", "up", "pair", "up", "empty"])
        item.set_style_payload(
            {
                "stroke_color": "#334455",
                "fill_color": "#f5e7a1",
                "label_color": "#1e40af",
                "arrow_up_color": "#000000",
                "arrow_down_color": "#d62828",
                "fill_visible": False,
                "box_base_visible": True,
            }
        )
        item.setRotation(18.0)
        item.setSelected(True)
        source._sync_selection_from_scene()

        assert source.apply_opacity_percent(65)

        payload = source._build_selection_payload()
        assert payload is not None
        assert len(payload["energy_diagrams"]) == 1
        assert payload["energy_diagrams"][0]["kind"] == "custom_level"
        assert payload["energy_diagrams"][0]["slot_count"] == 5
        assert payload["energy_diagrams"][0]["occupancies"] == ["pair", "up", "pair", "up", "empty"]
        assert payload["energy_diagrams"][0]["opacity"] == pytest.approx(0.65)

        target._last_scene_pos = QPointF(320.0, 260.0)
        target._paste_selection_payload(payload)

        assert len(target.energy_diagram_items) == 1
        pasted = target.energy_diagram_items[0]
        assert pasted.kind() == "custom_level"
        assert pasted.box_count() == 5
        assert pasted.label() == "4d"
        assert pasted.label_side() == "right"
        assert pasted.occupancies() == ("pair", "up", "pair", "up", "empty")
        assert pasted.style_payload() == {
            "stroke_color": "#334455",
            "fill_color": "#f5e7a1",
            "label_color": "#1e40af",
            "arrow_up_color": "#000000",
            "arrow_down_color": "#d62828",
            "fill_visible": False,
            "box_base_visible": True,
        }
        assert pasted.rotation() == pytest.approx(18.0)
        assert target.effective_item_opacity(pasted) == pytest.approx(0.65)
        assert pasted.isSelected() is True
    finally:
        source.close()
        target.close()


def test_energy_diagram_direct_keyboard_editing_uses_arrow_keys() -> None:
    canvas = ChemusonCanvas()
    try:
        canvas.resize(900, 700)
        canvas.show()
        QApplication.processEvents()

        item = canvas._insert_energy_diagram_item(QPointF(180.0, 140.0), kind="sublevel_p")

        assert item is not None
        second_slot = item.mapToScene(item._slot_regions()[1].center())
        QTest.mouseDClick(
            canvas.viewport(),
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
            canvas.mapFromScene(second_slot),
        )
        QApplication.processEvents()

        assert item.is_editing() is True
        assert item.occupancies() == ("empty", "empty", "empty")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Up)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "up", "empty")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Up)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "upup", "empty")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Down)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "pair", "empty")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Down)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "downdown", "empty")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Right)
        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Down)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "downdown", "down")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_A)
        QApplication.processEvents()
        assert item.occupancies() == ("empty", "downdown", "down")

        QTest.keyClick(canvas.viewport(), Qt.Key.Key_Escape)
        QApplication.processEvents()
        assert item.is_editing() is False
    finally:
        canvas.close()


def test_box_base_mode_switches_off_outline_immediately() -> None:
    canvas = ChemusonCanvas()
    try:
        item = canvas._insert_energy_diagram_item(QPointF(180.0, 140.0), kind="sublevel_d")

        assert item is not None
        item.setSelected(True)
        canvas._sync_selection_from_scene()

        assert item.effective_style()["box_stroke_visible"] is True
        assert item.effective_style()["box_base_visible"] is False

        canvas._set_selected_energy_diagram_box_base_visible(True)

        assert item.effective_style()["box_stroke_visible"] is False
        assert item.effective_style()["box_base_visible"] is True

        canvas._set_selected_energy_diagram_box_stroke_visible(True)

        assert item.effective_style()["box_stroke_visible"] is True
        assert item.effective_style()["box_base_visible"] is False
    finally:
        canvas.close()
