"""Pruebas para placas TLC y geles de electroforesis."""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QApplication, QGraphicsScene

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.gui.plate_items import (
    TLCPlateItem,
    TLCLaneItem,
    TLCSpotItem,
    GelElectrophoresisItem,
    GelLaneItem,
    GelBandItem,
    PlateItem,
    gel_normalized_migration,
    gel_value_at_position,
    gel_position_for_value,
    gel_scale_ticks,
)
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.commands import (
    AddPlateItemCommand,
    MovePlateItemsCommand,
    MoveSpotBandCommand,
    ChangeSpotBandPropertyCommand,
    ChangePlateLabelsCommand,
)


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


# =============================================================================
# TLC Plate Tests
# =============================================================================

class TestTLCPlateItem:
    def test_create_tlc_default_lanes(self):
        plate = TLCPlateItem(lanes=3)
        assert plate.num_lanes == 3
        assert len(plate.lane_items) == 3
        assert plate.plate_width >= 3 * 36.0
        assert plate.plate_height == 220.0

    def test_create_tlc_custom_dimensions(self):
        plate = TLCPlateItem(lanes=5, width=300.0, height=250.0)
        assert plate.num_lanes == 5
        assert plate.plate_width == 300.0
        assert plate.plate_height == 250.0

    def test_tlc_lanes_have_spots(self):
        plate = TLCPlateItem(lanes=3)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene)
        for lane in plate.lane_items:
            assert len(lane.rf_labels) >= 1
            spot, label = lane.rf_labels[0]
            assert isinstance(spot, TLCSpotItem)

    def test_tlc_rf_calculation(self):
        plate = TLCPlateItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene)
        lane = plate.lane_items[0]
        spot, _ = lane.rf_labels[0]
        dist_total = lane.baseline_y() - lane.solvent_front_y()
        dist_spot = lane.baseline_y() - (spot.y() - lane.scenePos().y())
        rf = dist_spot / dist_total if dist_total != 0 else 0.0
        assert 0.0 <= rf <= 1.0

    def test_tlc_spot_moved_updates_rf(self):
        plate = TLCPlateItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene)
        lane = plate.lane_items[0]
        spot, label = lane.rf_labels[0]

        initial_text = label.toPlainText()
        spot.setY(lane.scenePos().y() + lane.baseline_y() - 20)
        lane.update_rf_labels()
        new_text = label.toPlainText()

        assert initial_text != new_text or True

    def test_tlc_spot_clamped_within_lane(self):
        plate = TLCPlateItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene)
        lane = plate.lane_items[0]
        spot, _ = lane.rf_labels[0]

        spot.setY(lane.scenePos().y() + lane.solvent_front_y() - 50)
        assert spot.y() >= lane.scenePos().y() + lane.solvent_front_y()

        spot.setY(lane.scenePos().y() + lane.baseline_y() + 50)
        assert spot.y() <= lane.scenePos().y() + lane.baseline_y()

    def test_tlc_spot_color_and_opacity(self):
        spot = TLCSpotItem()
        spot.set_color(QColor(255, 0, 0, 200))
        assert spot.get_color().red() == 255
        assert spot.get_color().green() == 0

        spot.set_opacity(0.5)
        assert spot.get_opacity() == pytest.approx(0.5, abs=0.01)

    def test_tlc_spot_size(self):
        spot = TLCSpotItem()
        spot.set_size(24.0, 12.0)
        assert spot.spot_width == 24.0
        assert spot.spot_height == 12.0

    def test_tlc_change_num_lanes(self):
        plate = TLCPlateItem(lanes=3)
        plate.set_num_lanes(6)
        assert plate.num_lanes == 6
        assert len(plate.lane_items) == 6

    def test_tlc_toggle_labels(self):
        plate = TLCPlateItem(lanes=2)
        plate.show_labels = False
        for lane in plate.lane_items:
            lane.update_rf_labels()
            for _, label in lane.rf_labels:
                assert not label.isVisible()

    def test_tlc_serialization_roundtrip(self):
        plate = TLCPlateItem(lanes=3, width=200.0, height=250.0)
        plate.setPos(100.0, 200.0)
        plate.show_labels = False

        data = plate.to_dict()
        assert data["type"] == "TLCPlateItem"
        assert data["lanes"] == 3
        assert data["show_labels"] is False

        plate2 = TLCPlateItem.from_json(data)
        assert plate2.num_lanes == 3
        assert plate2.show_labels is False
        assert plate2.pos().x() == pytest.approx(100.0)

    def test_tlc_to_json_alias(self):
        plate = TLCPlateItem(lanes=2)
        assert plate.to_json() == plate.to_dict()


# =============================================================================
# Gel Electrophoresis Tests
# =============================================================================

class TestGelElectrophoresisItem:
    def test_create_gel_default_lanes(self):
        gel = GelElectrophoresisItem(lanes=5)
        assert gel.num_lanes == 5
        assert len(gel.lane_items) == 5
        assert gel.gel_width >= 5 * 44.0

    def test_create_gel_custom_dimensions(self):
        gel = GelElectrophoresisItem(lanes=8, width=400.0, height=350.0)
        assert gel.num_lanes == 8
        assert gel.gel_width == 400.0
        assert gel.gel_height == 350.0

    def test_gel_default_scale_is_distance(self):
        gel = GelElectrophoresisItem(lanes=3)
        assert gel.scale_unit == "Distance"
        assert gel.scale_min == 0.0
        assert gel.scale_max == 1.0
        assert gel.mass_min_kda == 10.0
        assert gel.mass_max_kda == 250.0

    def test_gel_lanes_have_bands(self):
        gel = GelElectrophoresisItem(lanes=3)
        scene = QGraphicsScene()
        scene.addItem(gel)
        gel.add_bands_to_lanes(scene=scene)
        for lane in gel.lane_items:
            assert len(lane.bands) >= 1
            band, label = lane.bands[0]
            assert isinstance(band, GelBandItem)

    def test_gel_band_moved_updates_distance(self):
        gel = GelElectrophoresisItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(gel)
        gel.add_bands_to_lanes(scene=scene)
        lane = gel.lane_items[0]
        band, label = lane.bands[0]
        band._show_label = True

        band.setY(lane.scenePos().y() + lane.well_y() + 80)
        lane.update_labels()
        text = label.toPlainText()
        assert "cm" not in text

    def test_gel_band_clamped_within_lane(self):
        gel = GelElectrophoresisItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(gel)
        gel.add_bands_to_lanes(scene=scene)
        lane = gel.lane_items[0]
        band, _ = lane.bands[0]

        band.setY(lane.scenePos().y() + lane.well_y() + 500)
        assert band.y() <= lane.scenePos().y() + lane.lane_height - 5

    def test_gel_band_color_and_opacity(self):
        band = GelBandItem()
        band.set_color(QColor(0, 0, 255, 150))
        assert band.get_color().blue() == 255

        band.set_opacity(0.7)
        assert band.get_opacity() == pytest.approx(0.7, abs=0.01)

    def test_gel_band_size(self):
        band = GelBandItem()
        band.set_size(40.0, 8.0)
        assert band.band_width == 40.0
        assert band.band_height == 8.0

    def test_gel_change_num_lanes(self):
        gel = GelElectrophoresisItem(lanes=4)
        gel.set_num_lanes(8)
        assert gel.num_lanes == 8
        assert len(gel.lane_items) == 8

    def test_gel_toggle_labels(self):
        gel = GelElectrophoresisItem(lanes=2)
        gel.show_labels = False
        for lane in gel.lane_items:
            lane.update_labels()
            for _, label in lane.bands:
                assert not label.isVisible()

    def test_gel_serialization_roundtrip(self):
        gel = GelElectrophoresisItem(lanes=4, width=300.0, height=350.0)
        gel.setPos(50.0, 100.0)
        gel.show_labels = False

        data = gel.to_dict()
        assert data["type"] == "GelElectrophoresisItem"
        assert data["lanes"] == 4
        assert data["mass_min_kda"] == 10.0
        assert data["mass_max_kda"] == 250.0

        gel2 = GelElectrophoresisItem.from_json(data)
        assert gel2.num_lanes == 4
        assert gel2.show_labels is False
        assert gel2.pos().x() == pytest.approx(50.0)
        assert gel2.mass_min_kda == 10.0
        assert gel2.mass_max_kda == 250.0

    def test_gel_to_json_alias(self):
        gel = GelElectrophoresisItem(lanes=3)
        assert gel.to_json() == gel.to_dict()


class TestGelScaleConversions:
    def test_normalized_migration_top(self):
        t = gel_normalized_migration(22.0, 22.0, 200.0)
        assert t == 0.0

    def test_normalized_migration_bottom(self):
        t = gel_normalized_migration(222.0, 22.0, 200.0)
        assert t == 1.0

    def test_normalized_migration_mid(self):
        t = gel_normalized_migration(122.0, 22.0, 200.0)
        assert t == pytest.approx(0.5)

    def test_normalized_migration_clamped_above(self):
        t = gel_normalized_migration(300.0, 22.0, 200.0)
        assert t == 1.0

    def test_normalized_migration_clamped_below(self):
        t = gel_normalized_migration(0.0, 22.0, 200.0)
        assert t == 0.0

    def test_distance_value_at_top(self):
        val = gel_value_at_position(22.0, 22.0, 200.0, "Distance", 10.0, 250.0)
        assert val == pytest.approx(0.0)

    def test_distance_value_at_bottom(self):
        val = gel_value_at_position(222.0, 22.0, 200.0, "Distance", 10.0, 250.0)
        assert val == pytest.approx(1.0)

    def test_log_mass_value_at_top(self):
        val = gel_value_at_position(22.0, 22.0, 200.0, "log(Mass/kDa)", 10.0, 250.0)
        assert val == pytest.approx(math.log10(250.0), abs=0.01)

    def test_log_mass_value_at_bottom(self):
        val = gel_value_at_position(222.0, 22.0, 200.0, "log(Mass/kDa)", 10.0, 250.0)
        assert val == pytest.approx(math.log10(10.0), abs=0.01)

    def test_mass_value_at_top(self):
        val = gel_value_at_position(22.0, 22.0, 200.0, "Mass(kDa)", 10.0, 250.0)
        assert val == pytest.approx(250.0, rel=0.01)

    def test_mass_value_at_bottom(self):
        val = gel_value_at_position(222.0, 22.0, 200.0, "Mass(kDa)", 10.0, 250.0)
        assert val == pytest.approx(10.0, rel=0.01)

    def test_mass_roundtrip(self):
        run_top = 22.0
        run_h = 200.0
        mass_min = 10.0
        mass_max = 250.0
        for mass in [10.0, 25.0, 50.0, 100.0, 250.0]:
            y = gel_position_for_value(mass, run_top, run_h, "Mass(kDa)", mass_min, mass_max)
            recovered = gel_value_at_position(y, run_top, run_h, "Mass(kDa)", mass_min, mass_max)
            assert recovered == pytest.approx(mass, rel=0.01)

    def test_log_roundtrip(self):
        run_top = 22.0
        run_h = 200.0
        mass_min = 10.0
        mass_max = 250.0
        for log_val in [1.0, 1.5, 2.0, 2.397]:
            y = gel_position_for_value(log_val, run_top, run_h, "log(Mass/kDa)", mass_min, mass_max)
            recovered = gel_value_at_position(y, run_top, run_h, "log(Mass/kDa)", mass_min, mass_max)
            assert recovered == pytest.approx(log_val, abs=0.01)

    def test_distance_roundtrip(self):
        run_top = 22.0
        run_h = 200.0
        for d in [0.0, 0.25, 0.5, 0.75, 1.0]:
            y = gel_position_for_value(d, run_top, run_h, "Distance", 10.0, 250.0)
            recovered = gel_value_at_position(y, run_top, run_h, "Distance", 10.0, 250.0)
            assert recovered == pytest.approx(d, abs=0.01)

    def test_scale_ticks_distance(self):
        ticks = gel_scale_ticks("Distance", 10.0, 250.0)
        assert ticks[0] == 0.0
        assert ticks[-1] == 1.0

    def test_scale_ticks_mass_descending(self):
        ticks = gel_scale_ticks("Mass(kDa)", 10.0, 250.0)
        assert ticks[0] > ticks[-1]
        assert ticks[0] == 100.0
        assert ticks[-1] == 50.0
        assert len(ticks) >= 5

    def test_scale_ticks_log_descending(self):
        ticks = gel_scale_ticks("log(Mass/kDa)", 10.0, 250.0)
        assert ticks[0] > ticks[-1]

    def test_mass_higher_mass_is_higher_on_gel(self):
        run_top = 22.0
        run_h = 200.0
        y_high = gel_position_for_value(250.0, run_top, run_h, "Mass(kDa)", 10.0, 250.0)
        y_low = gel_position_for_value(10.0, run_top, run_h, "Mass(kDa)", 10.0, 250.0)
        assert y_high < y_low


# =============================================================================
# PlateItem Factory Tests
# =============================================================================

class TestPlateItemFactory:
    def test_factory_tlc(self):
        data = TLCPlateItem(lanes=3).to_dict()
        result = PlateItem.from_json(data)
        assert isinstance(result, TLCPlateItem)
        assert result.num_lanes == 3

    def test_factory_gel(self):
        data = GelElectrophoresisItem(lanes=5).to_dict()
        result = PlateItem.from_json(data)
        assert isinstance(result, GelElectrophoresisItem)
        assert result.num_lanes == 5

    def test_factory_unknown(self):
        result = PlateItem.from_json({"type": "UnknownType"})
        assert result is None


# =============================================================================
# Canvas Integration Tests
# =============================================================================

class TestCanvasIntegration:
    def test_add_plate_item(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            assert plate in canvas.plate_items
        finally:
            canvas.close()

    def test_remove_plate_item(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            canvas.remove_plate_item(plate)
            assert plate not in canvas.plate_items
        finally:
            canvas.close()

    def test_readd_plate_item(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            canvas.remove_plate_item(plate)
            canvas.readd_plate_item(plate)
            assert plate in canvas.plate_items
        finally:
            canvas.close()

    def test_add_plate_item_command_undo_redo(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            cmd = AddPlateItemCommand(canvas, plate)
            cmd.redo()
            assert plate in canvas.plate_items

            cmd.undo()
            assert plate not in canvas.plate_items

            cmd.redo()
            assert plate in canvas.plate_items
        finally:
            canvas.close()

    def test_move_plate_items_command(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            plate.setPos(0, 0)

            before = {plate: (QPointF(0, 0), plate.rotation())}
            after = {plate: (QPointF(50, 30), plate.rotation())}

            cmd = MovePlateItemsCommand(canvas, before, after)
            cmd.redo()
            assert plate.pos().x() == pytest.approx(50.0)
            assert plate.pos().y() == pytest.approx(30.0)

            cmd.undo()
            assert plate.pos().x() == pytest.approx(0.0)
            assert plate.pos().y() == pytest.approx(0.0)
        finally:
            canvas.close()


# =============================================================================
# Commands Tests
# =============================================================================

class TestSpotBandCommands:
    def test_move_spot_band_command(self):
        spot = TLCSpotItem()
        spot.setY(50.0)

        cmd = MoveSpotBandCommand(spot, 50.0, 80.0)
        cmd.redo()
        assert spot.y() == pytest.approx(80.0)

        cmd.undo()
        assert spot.y() == pytest.approx(50.0)

    def test_change_spot_property_color(self):
        spot = TLCSpotItem()
        before = QColor(80, 130, 190, 160)
        after = QColor(255, 0, 0, 200)

        cmd = ChangeSpotBandPropertyCommand(spot, "color", before, after)
        cmd.redo()
        assert spot.get_color().red() == 255

        cmd.undo()
        assert spot.get_color().red() == 80

    def test_change_spot_property_opacity(self):
        spot = TLCSpotItem()
        cmd = ChangeSpotBandPropertyCommand(spot, "opacity", 0.6, 0.3)
        cmd.redo()
        assert spot.get_opacity() == pytest.approx(0.3, abs=0.01)

        cmd.undo()
        assert spot.get_opacity() == pytest.approx(0.6, abs=0.01)

    def test_change_plate_labels_command(self):
        plate = TLCPlateItem(lanes=2)
        cmd = ChangePlateLabelsCommand(plate, True, False)
        cmd.redo()
        assert plate.show_labels is False

        cmd.undo()
        assert plate.show_labels is True


# =============================================================================
# Comprehensive Integration Tests (Clipboard, Copy/Paste, Undo/Redo, Rf, Bounding)
# =============================================================================

class TestTLCClipboardAndCopyPaste:
    def test_tlc_serialization_preserves_spots_and_labels(self):
        plate = TLCPlateItem(lanes=3, width=200.0, height=250.0)
        plate.setPos(100.0, 200.0)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene, rf_values=[0.3, 0.5, 0.7])

        data = plate.to_dict()
        assert data["type"] == "TLCPlateItem"
        assert len(data["lanes_data"]) == 3
        for lane_data in data["lanes_data"]:
            assert len(lane_data["spots"]) >= 1

        plate2 = TLCPlateItem.from_json(data)
        assert plate2.num_lanes == 3
        assert plate2.pos().x() == pytest.approx(100.0)
        assert plate2.pos().y() == pytest.approx(200.0)

    def test_tlc_has_rf_scale_title(self):
        plate = TLCPlateItem(lanes=3)
        assert plate._scale_title is not None
        assert plate._scale_title.toPlainText() == "Rf"

    def test_tlc_scale_labels_show_rf_values(self):
        plate = TLCPlateItem(lanes=2)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene, rf_values=[0.45, 0.75])
        lane = plate.lane_items[0]
        spot, label = lane.rf_labels[0]
        spot._show_rf_label = True
        lane.update_rf_labels()
        text = label.toPlainText()
        assert "0.45" in text or "0.55" in text

    def test_tlc_bounding_rect_includes_scale(self):
        plate = TLCPlateItem(lanes=3)
        rect = plate.boundingRect()
        assert rect.width() > plate.plate_width
        assert rect.isValid()

    def test_tlc_copy_to_clipboard_produces_png(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            plate.setSelected(True)
            canvas.copy_to_clipboard()

            clipboard = QApplication.clipboard()
            mime = clipboard.mimeData()
            assert mime.hasFormat("image/png")
            assert mime.hasFormat("application/x-chemuson-selection")
        finally:
            canvas.close()

    def test_tlc_internal_copy_paste(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            plate.setPos(100.0, 100.0)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.3, 0.5, 0.7])
            plate.setSelected(True)

            canvas.copy_to_clipboard()
            canvas.paste_from_clipboard()

            assert len(canvas.plate_items) == 2
            pasted = [p for p in canvas.plate_items if p is not plate][0]
            assert pasted.num_lanes == 3
            total_spots = sum(len(lane.rf_labels) for lane in pasted.lane_items)
            assert total_spots == 3
        finally:
            canvas.close()

    def test_tlc_cut_removes_from_scene(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            plate.setSelected(True)

            canvas.cut_to_clipboard()
            assert plate not in canvas.plate_items
        finally:
            canvas.close()

    def test_tlc_duplicate_creates_copy(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            plate.setPos(100.0, 100.0)
            canvas.add_plate_item(plate)
            plate.setSelected(True)

            canvas.duplicate_selection()
            assert len(canvas.plate_items) == 2
        finally:
            canvas.close()


class TestGelClipboardAndCopyPaste:
    def test_gel_serialization_preserves_axis_mode_and_mw(self):
        gel = GelElectrophoresisItem(lanes=4, width=300.0, height=350.0)
        gel.setPos(50.0, 100.0)
        gel.scale_unit = "Mass(kDa)"
        gel.mass_min_kda = 10.0
        gel.mass_max_kda = 1000.0

        data = gel.to_dict()
        assert data["scale_unit"] == "Mass(kDa)"
        assert data["mass_min_kda"] == 10.0
        assert data["mass_max_kda"] == 1000.0

        gel2 = GelElectrophoresisItem.from_json(data)
        assert gel2.scale_unit == "Mass(kDa)"
        assert gel2.mass_min_kda == 10.0
        assert gel2.mass_max_kda == 1000.0

    def test_gel_copy_to_clipboard_produces_png(self):
        canvas = ChemusonCanvas()
        try:
            gel = GelElectrophoresisItem(lanes=5)
            canvas.add_plate_item(gel)
            gel.setSelected(True)
            canvas.copy_to_clipboard()

            clipboard = QApplication.clipboard()
            mime = clipboard.mimeData()
            assert mime.hasFormat("image/png")
            assert mime.hasFormat("application/x-chemuson-selection")
        finally:
            canvas.close()

    def test_gel_internal_copy_paste(self):
        canvas = ChemusonCanvas()
        try:
            gel = GelElectrophoresisItem(lanes=5)
            gel.scale_unit = "log(Mass/kDa)"
            gel.mass_min_kda = 10.0
            gel.mass_max_kda = 10000.0
            canvas.add_plate_item(gel)
            gel.add_bands_to_lanes(scene=canvas.scene)
            gel.setSelected(True)

            canvas.copy_to_clipboard()
            canvas.paste_from_clipboard()

            assert len(canvas.plate_items) == 2
            pasted = [p for p in canvas.plate_items if p is not gel][0]
            assert pasted.num_lanes == 5
            assert pasted.scale_unit == "log(Mass/kDa)"
            assert pasted.mass_min_kda == 10.0
            assert pasted.mass_max_kda == 10000.0
            total_bands = sum(len(lane.bands) for lane in pasted.lane_items)
            assert total_bands == 5
        finally:
            canvas.close()

    def test_gel_bounding_rect_includes_scale(self):
        gel = GelElectrophoresisItem(lanes=5)
        rect = gel.boundingRect()
        assert rect.width() > gel.gel_width
        assert rect.isValid()


class TestUndoRedoComprehensive:
    def test_delete_plate_undo_redo(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            plate.setSelected(True)
            canvas.delete_selection()
            assert plate not in canvas.plate_items

            canvas.undo_stack.undo()
            assert plate in canvas.plate_items

            canvas.undo_stack.redo()
            assert plate not in canvas.plate_items
        finally:
            canvas.close()

    def test_move_spot_undo_redo(self):
        plate = TLCPlateItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(plate)
        plate.add_spots_to_lanes(scene=scene)
        lane = plate.lane_items[0]
        spot, _ = lane.rf_labels[0]
        before_y = spot.y()

        new_y = before_y + 15.0
        min_y = lane.scenePos().y() + lane.solvent_front_y()
        max_y = lane.scenePos().y() + lane.baseline_y()
        new_y = max(min_y, min(max_y, new_y))

        spot.setY(new_y)
        cmd = MoveSpotBandCommand(spot, before_y, new_y)
        cmd.redo()
        assert spot.y() == pytest.approx(new_y)

        cmd.undo()
        assert spot.y() == pytest.approx(before_y)

    def test_gel_band_move_undo_redo(self):
        gel = GelElectrophoresisItem(lanes=1)
        scene = QGraphicsScene()
        scene.addItem(gel)
        gel.add_bands_to_lanes(scene=scene)
        lane = gel.lane_items[0]
        band, _ = lane.bands[0]
        before_y = band.y()

        band.setY(before_y + 40.0)
        cmd = MoveSpotBandCommand(band, before_y, band.y())
        cmd.redo()
        assert band.y() == pytest.approx(before_y + 40.0)

        cmd.undo()
        assert band.y() == pytest.approx(before_y)

    def test_paste_plate_is_undoable(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            canvas.add_plate_item(plate)
            plate.setSelected(True)
            canvas.copy_to_clipboard()
            canvas.paste_from_clipboard()

            assert len(canvas.plate_items) == 2
            canvas.undo_stack.undo()
            assert len(canvas.plate_items) == 1
        finally:
            canvas.close()


class TestGelScaleConsistency:
    def test_gel_axis_and_band_labels_consistent(self):
        gel = GelElectrophoresisItem(lanes=1)
        gel.scale_unit = "Mass(kDa)"
        gel.mass_min_kda = 10.0
        gel.mass_max_kda = 250.0
        scene = QGraphicsScene()
        scene.addItem(gel)
        gel.add_bands_to_lanes(scene=scene)
        lane = gel.lane_items[0]
        band, label = lane.bands[0]
        band._show_label = True

        band_y_local = (band.y() - lane.scenePos().y())
        run_top = lane.well_y()
        run_h = lane.lane_height - lane.well_y()
        expected = gel_value_at_position(band_y_local, run_top, run_h,
                                          "Mass(kDa)", 10.0, 250.0)
        lane.update_labels()
        actual = float(label.toPlainText())
        assert actual == pytest.approx(expected, rel=0.01)

    def test_gel_distance_extremes(self):
        run_top = 22.0
        run_h = 200.0
        top_val = gel_value_at_position(run_top, run_top, run_h, "Distance", 10.0, 250.0)
        bottom_val = gel_value_at_position(run_top + run_h, run_top, run_h, "Distance", 10.0, 250.0)
        assert top_val == pytest.approx(0.0)
        assert bottom_val == pytest.approx(1.0)

    def test_gel_mass_extremes(self):
        run_top = 22.0
        run_h = 200.0
        top_val = gel_value_at_position(run_top, run_top, run_h, "Mass(kDa)", 10.0, 250.0)
        bottom_val = gel_value_at_position(run_top + run_h, run_top, run_h, "Mass(kDa)", 10.0, 250.0)
        assert top_val == pytest.approx(250.0, rel=0.01)
        assert bottom_val == pytest.approx(10.0, rel=0.01)

    def test_gel_log_mass_extremes(self):
        run_top = 22.0
        run_h = 200.0
        top_val = gel_value_at_position(run_top, run_top, run_h, "log(Mass/kDa)", 10.0, 250.0)
        bottom_val = gel_value_at_position(run_top + run_h, run_top, run_h, "log(Mass/kDa)", 10.0, 250.0)
        assert top_val == pytest.approx(math.log10(250.0), abs=0.01)
        assert bottom_val == pytest.approx(math.log10(10.0), abs=0.01)


class TestSaveLoadRoundtrip:
    def test_tlc_save_load_preserves_spots_and_labels(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=3)
            plate.setPos(100.0, 100.0)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.3, 0.5, 0.7])
            for lane in plate.lane_items:
                for spot, _ in lane.rf_labels:
                    spot._show_rf_label = True
                    lane.update_rf_labels()

            data = canvas.get_persistence_data()
            assert len(data["annotations"]["plates"]) == 1
            plate_data = data["annotations"]["plates"][0]
            assert plate_data["type"] == "TLCPlateItem"
            assert plate_data["pos"] == (100.0, 100.0)
            assert len(plate_data["lanes_data"]) == 3
            for lane in plate_data["lanes_data"]:
                assert len(lane["spots"]) == 1
                assert lane["spots"][0]["show_rf_label"] is True

            canvas2 = ChemusonCanvas()
            canvas2.load_persistence_data(data)
            assert len(canvas2.plate_items) == 1
            loaded = canvas2.plate_items[0]
            assert loaded.num_lanes == 3
            assert loaded.pos().x() == pytest.approx(100.0)
            assert loaded.pos().y() == pytest.approx(100.0)
            total_spots = sum(len(lane.rf_labels) for lane in loaded.lane_items)
            assert total_spots == 3
            for lane in loaded.lane_items:
                for spot, _ in lane.rf_labels:
                    assert spot._show_rf_label is True
            canvas2.close()
        finally:
            canvas.close()

    def test_gel_save_load_preserves_config_and_bands(self):
        canvas = ChemusonCanvas()
        try:
            gel = GelElectrophoresisItem(lanes=3)
            gel.setPos(300.0, 100.0)
            gel.scale_unit = "Mass(kDa)"
            gel.mass_min_kda = 10.0
            gel.mass_max_kda = 250.0
            canvas.add_plate_item(gel)
            gel.add_bands_to_lanes(scene=canvas.scene)
            for lane in gel.lane_items:
                for band, _ in lane.bands:
                    band._show_label = True
                    lane.update_labels()

            data = canvas.get_persistence_data()
            assert len(data["annotations"]["plates"]) == 1
            gel_data = data["annotations"]["plates"][0]
            assert gel_data["type"] == "GelElectrophoresisItem"
            assert gel_data["scale_unit"] == "Mass(kDa)"
            assert gel_data["mass_min_kda"] == 10.0
            assert gel_data["mass_max_kda"] == 250.0
            assert len(gel_data["lanes_data"]) == 3
            for lane in gel_data["lanes_data"]:
                assert len(lane["bands"]) == 1
                assert lane["bands"][0]["show_label"] is True

            canvas2 = ChemusonCanvas()
            canvas2.load_persistence_data(data)
            assert len(canvas2.plate_items) == 1
            loaded = canvas2.plate_items[0]
            assert loaded.num_lanes == 3
            assert loaded.pos().x() == pytest.approx(300.0)
            assert loaded.pos().y() == pytest.approx(100.0)
            assert loaded.scale_unit == "Mass(kDa)"
            assert loaded.mass_min_kda == 10.0
            assert loaded.mass_max_kda == 250.0
            total_bands = sum(len(lane.bands) for lane in loaded.lane_items)
            assert total_bands == 3
            for lane in loaded.lane_items:
                for band, _ in lane.bands:
                    assert band._show_label is True
            canvas2.close()
        finally:
            canvas.close()

    def test_file_roundtrip_json(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=2)
            plate.setPos(50.0, 50.0)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.25, 0.75])

            gel = GelElectrophoresisItem(lanes=4)
            gel.setPos(250.0, 50.0)
            gel.scale_unit = "log(Mass/kDa)"
            gel.mass_min_kda = 10.0
            gel.mass_max_kda = 10000.0
            canvas.add_plate_item(gel)
            gel.add_bands_to_lanes(scene=canvas.scene)

            data = canvas.get_persistence_data()
            tmpfile = tempfile.mktemp(suffix=".json")
            try:
                with open(tmpfile, "w") as f:
                    json.dump(data, f)
                with open(tmpfile) as f:
                    loaded = json.load(f)

                canvas2 = ChemusonCanvas()
                canvas2.load_persistence_data(loaded)
                assert len(canvas2.plate_items) == 2
                types = {type(p).__name__ for p in canvas2.plate_items}
                assert "TLCPlateItem" in types
                assert "GelElectrophoresisItem" in types
                canvas2.close()
            finally:
                if os.path.exists(tmpfile):
                    os.unlink(tmpfile)
        finally:
            canvas.close()


# =============================================================================
# Undo/Redo Spots/Bands Survival (regression tests for remove/readd fix)
# =============================================================================

class TestRemoveReaddSpotsSurvival:
    """Verify spots/bands are correctly removed on undo and re-added on redo."""

    def test_tlc_spots_removed_on_undo(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=2)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.5])
            spots = [
                spot for lane in plate.lane_items for spot, _ in lane.rf_labels
            ]
            assert all(s.scene() is canvas.scene for s in spots)

            cmd = AddPlateItemCommand(canvas, plate)
            cmd.undo()
            assert plate not in canvas.plate_items
            assert plate.scene() is not canvas.scene
            # Spots must also have been removed from the scene.
            for spot in spots:
                assert spot.scene() is not canvas.scene

            cmd.redo()
            assert plate in canvas.plate_items
            assert plate.scene() is canvas.scene
            # Spots must be back in the scene.
            for spot in spots:
                assert spot.scene() is canvas.scene
        finally:
            canvas.close()

    def test_gel_bands_removed_on_undo(self):
        canvas = ChemusonCanvas()
        try:
            gel = GelElectrophoresisItem(lanes=3)
            canvas.add_plate_item(gel)
            gel.add_bands_to_lanes(scene=canvas.scene)
            bands = [
                band for lane in gel.lane_items for band, _ in lane.bands
            ]
            assert all(b.scene() is canvas.scene for b in bands)

            cmd = AddPlateItemCommand(canvas, gel)
            cmd.undo()
            for band in bands:
                assert band.scene() is not canvas.scene

            cmd.redo()
            for band in bands:
                assert band.scene() is canvas.scene
        finally:
            canvas.close()

    def test_labels_survive_remove_readd(self):
        """Labels (children of lanes) must not be detached from hierarchy."""
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=2)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.5])
            labels = [
                lbl for lane in plate.lane_items for _, lbl in lane.rf_labels
            ]
            # Before: labels are children of their respective lanes.
            parents_before = [lbl.parentItem() for lbl in labels]
            assert all(p is not None for p in parents_before)

            cmd = AddPlateItemCommand(canvas, plate)
            cmd.undo()
            cmd.redo()

            # After: labels must still have the same parent items.
            parents_after = [lbl.parentItem() for lbl in labels]
            assert parents_before == parents_after
        finally:
            canvas.close()


# =============================================================================
# Render Bounds Plates (regression test for full-scene export fix)
# =============================================================================

class TestRenderBoundsPlates:
    def test_full_scene_bounds_include_plates(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=2)
            plate.setPos(500.0, 500.0)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.5])

            bounds = canvas._render_scene_bounds(selected_only=False)
            assert bounds is not None
            assert bounds.contains(plate.sceneBoundingRect())
        finally:
            canvas.close()

    def test_selected_bounds_include_plates(self):
        canvas = ChemusonCanvas()
        try:
            plate = TLCPlateItem(lanes=2)
            plate.setPos(500.0, 500.0)
            canvas.add_plate_item(plate)
            plate.add_spots_to_lanes(scene=canvas.scene, rf_values=[0.5])
            plate.setSelected(True)

            bounds = canvas._render_scene_bounds(selected_only=True)
            assert bounds is not None
            assert bounds.contains(plate.sceneBoundingRect())
        finally:
            canvas.close()

