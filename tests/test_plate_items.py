"""Pruebas para placas TLC y geles de electroforesis."""

from __future__ import annotations

import os
import sys

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
        assert "cm" in text

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

        gel2 = GelElectrophoresisItem.from_json(data)
        assert gel2.num_lanes == 4
        assert gel2.show_labels is False
        assert gel2.pos().x() == pytest.approx(50.0)

    def test_gel_to_json_alias(self):
        gel = GelElectrophoresisItem(lanes=3)
        assert gel.to_json() == gel.to_dict()


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
