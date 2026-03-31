from __future__ import annotations
import math
from PyQt6.QtWidgets import (
    QGraphicsObject,
    QGraphicsItem,
    QGraphicsRectItem,
    QGraphicsEllipseItem,
    QGraphicsLineItem,
    QGraphicsTextItem,
    QGraphicsSceneMouseEvent,
    QGraphicsPathItem
)
from PyQt6.QtGui import QColor, QPen, QBrush, QFont, QPainter, QPainterPath
from PyQt6.QtCore import Qt, QRectF, QPointF, pyqtSignal

class TLCSpotItem(QGraphicsObject):
    moved = pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        
        self.spot_width = 16.0
        self.spot_height = 8.0
        self._color = QColor(100, 150, 200, 180) # Default color, mildly opaque
        
        self.spot_rect = QGraphicsEllipseItem(
            -self.spot_width/2, -self.spot_height/2,
            self.spot_width, self.spot_height, self
        )
        self.spot_rect.setPen(QPen(Qt.PenStyle.NoPen))
        self.spot_rect.setBrush(QBrush(self._color))
        self._update_appearance()

    def set_color(self, color: QColor):
        self._color = color
        self._update_appearance()
        
    def get_color(self) -> QColor:
        return self._color
        
    def set_opacity(self, opacity: float):
        alpha = int(opacity * 255)
        self._color.setAlpha(alpha)
        self._update_appearance()

    def _update_appearance(self):
        self.spot_rect.setBrush(QBrush(self._color))

    def boundingRect(self) -> QRectF:
        return self.spot_rect.boundingRect()
        
    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(Qt.GlobalColor.blue, 1, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect())

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange and self.parentItem():
            parent: TLCLaneItem = self.parentItem()
            new_pos = value
            new_pos.setX(0)
            y = new_pos.y()
            min_y = parent.solvent_front_y()
            max_y = parent.baseline_y()
            if y < min_y:
                y = min_y
            elif y > max_y:
                y = max_y
            new_pos.setY(y)
            return new_pos
        elif change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.moved.emit()
            if hasattr(self.parentItem(), 'update_rf_labels'):
                self.parentItem().update_rf_labels()
        return super().itemChange(change, value)


class TLCLaneItem(QGraphicsObject):
    def __init__(self, parent=None, index=0, width=40.0, height=200.0):
        super().__init__(parent)
        self.index = index
        self.lane_width = width
        self.lane_height = height
        self._baseline_y = height * 0.85
        self._solvent_front_y = height * 0.15
        
        self.rf_labels = []

    def set_boundaries(self, baseline_y: float, solvent_front_y: float):
        self._baseline_y = baseline_y
        self._solvent_front_y = solvent_front_y
        self.update_rf_labels()

    def baseline_y(self): return self._baseline_y
    def solvent_front_y(self): return self._solvent_front_y

    def add_spot(self, rf: float = 0.5):
        spot = TLCSpotItem(self)
        y = self._baseline_y - rf * (self._baseline_y - self._solvent_front_y)
        spot.setPos(0, y)
        
        label = QGraphicsTextItem("", self)
        font = QFont("Arial", 8)
        label.setFont(font)
        label.setDefaultTextColor(QColor("black"))
        self.rf_labels.append((spot, label))
        self.update_rf_labels()
        return spot

    def update_rf_labels(self):
        plate = self.parentItem()
        show_labels = True
        if hasattr(plate, 'show_labels'):
            show_labels = plate.show_labels

        dist_total = self._baseline_y - self._solvent_front_y
        for spot, label in self.rf_labels:
            if not show_labels:
                label.setVisible(False)
                continue
            label.setVisible(True)
            dist_spot = self._baseline_y - spot.y()
            rf = dist_spot / dist_total if dist_total != 0 else 0
            rf = max(0.0, min(1.0, rf))
            label.setPlainText(f"Rf={rf:.2f}")
            label.setPos(10, spot.y() - label.boundingRect().height()/2)

    def boundingRect(self) -> QRectF:
        return QRectF(-self.lane_width/2, 0, self.lane_width, self.lane_height)

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(100, 100, 255, 100), 1, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect())


class TLCPlateItem(QGraphicsObject):
    def __init__(self, lanes: int = 3, width: float = 150.0, height: float = 200.0, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        
        self.num_lanes = lanes
        self.plate_width = max(width, lanes * 30.0)
        self.plate_height = height
        
        self.baseline_y = self.plate_height * 0.85
        self.solvent_front_y = self.plate_height * 0.15
        self.show_labels = True

        self._setup_plate()

    def _setup_plate(self):
        self.outline = QGraphicsRectItem(0, 0, self.plate_width, self.plate_height, self)
        self.outline.setPen(QPen(Qt.GlobalColor.black, 2))
        self.outline.setBrush(QBrush(QColor(255, 255, 255, 230)))
        
        pen_dash = QPen(Qt.GlobalColor.gray, 1, Qt.PenStyle.DashLine)
        self.baseline_line = QGraphicsLineItem(10, self.baseline_y, self.plate_width-10, self.baseline_y, self)
        self.baseline_line.setPen(pen_dash)
        
        self.solvent_line = QGraphicsLineItem(10, self.solvent_front_y, self.plate_width-10, self.solvent_front_y, self)
        self.solvent_line.setPen(pen_dash)

        self.lane_items = []
        lane_w = self.plate_width / self.num_lanes
        for i in range(self.num_lanes):
            lane = TLCLaneItem(self, index=i, width=lane_w, height=self.plate_height)
            lane.setPos(lane_w * i + lane_w/2, 0)
            lane.set_boundaries(self.baseline_y, self.solvent_front_y)
            lane.add_spot(rf=0.2 + 0.1*(i%3))
            self.lane_items.append(lane)

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, self.plate_width, self.plate_height)
        
    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(Qt.GlobalColor.blue, 2, Qt.PenStyle.DotLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect().adjusted(-2, -2, 2, 2))

    def to_dict(self):
        lanes_data = []
        for lane in self.lane_items:
            spots = []
            for spot, _ in lane.rf_labels:
                spots.append({
                    "y": spot.y(),
                    "color": spot.get_color().name(QColor.NameFormat.HexArgb)
                })
            lanes_data.append(spots)
            
        return {
            "type": "TLCPlateItem",
            "pos": (self.pos().x(), self.pos().y()),
            "lanes": self.num_lanes,
            "width": self.plate_width,
            "height": self.plate_height,
            "show_labels": self.show_labels,
            "lanes_data": lanes_data
        }

    def load_dict(self, data):
        self.setPos(data.get("pos", (0,0))[0], data.get("pos", (0,0))[1])
        self.show_labels = data.get("show_labels", True)
        lanes_data = data.get("lanes_data", [])
        for i, lane_data in enumerate(lanes_data):
            if i < len(self.lane_items):
                lane = self.lane_items[i]
                for spot, lbl in lane.rf_labels:
                    self.scene().removeItem(spot)
                    self.scene().removeItem(lbl)
                lane.rf_labels.clear()
                
                for sdata in lane_data:
                    spot = TLCSpotItem(lane)
                    spot.setPos(0, sdata.get("y", lane._baseline_y))
                    if "color" in sdata:
                        spot.set_color(QColor(sdata["color"]))
                    
                    lbl = QGraphicsTextItem("", lane)
                    lbl.setFont(QFont("Arial", 8))
                    lane.rf_labels.append((spot, lbl))
                lane.update_rf_labels()

    def to_json(self) -> dict:
        return self.to_dict()

    @classmethod
    def from_json(cls, data: dict) -> "PlateItem":
        item = TLCPlateItem(lanes=data.get("lanes", 3), width=data.get("width", 150.0), height=data.get("height", 200.0))
        item.load_dict(data)
        return item


class GelBandItem(QGraphicsObject):
    moved = pyqtSignal()
    
    def __init__(self, parent=None, width=30.0):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        
        self.band_width = width
        self.band_height = 6.0
        self._color = QColor(50, 50, 50, 200) # Dark grey, somewhat opaque
        
        self.band_rect = QGraphicsRectItem(
            -self.band_width/2, -self.band_height/2,
            self.band_width, self.band_height, self
        )
        self.band_rect.setPen(QPen(Qt.PenStyle.NoPen))
        self.band_rect.setBrush(QBrush(self._color))

    def set_color(self, color: QColor):
        self._color = color
        self.band_rect.setBrush(QBrush(self._color))

    def get_color(self) -> QColor:
        return self._color

    def set_opacity(self, opacity: float):
        alpha = int(opacity * 255)
        self._color.setAlpha(alpha)
        self.band_rect.setBrush(QBrush(self._color))

    def boundingRect(self) -> QRectF:
        return self.band_rect.boundingRect()
        
    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(Qt.GlobalColor.blue, 1, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect())

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange and self.parentItem():
            parent: GelLaneItem = self.parentItem()
            new_pos = value
            new_pos.setX(0)
            y = new_pos.y()
            min_y = parent.well_y()
            max_y = parent.lane_height
            if y < min_y:
                y = min_y
            elif y > max_y:
                y = max_y
            new_pos.setY(y)
            return new_pos
        elif change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.moved.emit()
            if hasattr(self.parentItem(), 'update_labels'):
                self.parentItem().update_labels()
        return super().itemChange(change, value)


class GelLaneItem(QGraphicsObject):
    def __init__(self, parent=None, index=0, width=40.0, height=250.0):
        super().__init__(parent)
        self.index = index
        self.lane_width = width
        self.lane_height = height
        self._well_y = 20.0
        
        self.bands = []

    def well_y(self): return self._well_y

    def add_band(self, dist: float = 50.0):
        band = GelBandItem(self, width=self.lane_width*0.8)
        band.setPos(0, self._well_y + dist)
        
        label = QGraphicsTextItem("", self)
        font = QFont("Arial", 8)
        label.setFont(font)
        self.bands.append((band, label))
        self.update_labels()
        return band

    def update_labels(self):
        gel = self.parentItem()
        show_labels = True
        if hasattr(gel, 'show_labels'):
            show_labels = gel.show_labels

        for band, label in self.bands:
            if not show_labels:
                label.setVisible(False)
                continue
            label.setVisible(True)
            dist = band.y() - self._well_y
            label.setPlainText(f"{dist*0.05:.1f} cm")
            label.setPos(self.lane_width/2 + 2, band.y() - label.boundingRect().height()/2)

    def boundingRect(self) -> QRectF:
        return QRectF(-self.lane_width/2, 0, self.lane_width, self.lane_height)

    def paint(self, painter, option, widget=None):
        pass


class GelElectrophoresisItem(QGraphicsObject):
    def __init__(self, lanes: int = 5, width: float = 250.0, height: float = 300.0, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        
        self.num_lanes = lanes
        self.gel_width = max(width, lanes * 40.0)
        self.gel_height = height
        self.show_labels = True
        
        self._setup_gel()

    def _setup_gel(self):
        self.outline = QGraphicsRectItem(0, 0, self.gel_width, self.gel_height, self)
        self.outline.setPen(QPen(Qt.GlobalColor.darkGray, 1))
        self.outline.setBrush(QBrush(QColor(240, 240, 245, 200)))
        
        self.lane_items = []
        lane_w = self.gel_width / self.num_lanes
        well_height = 12.0
        well_width = lane_w * 0.7
        
        path = QPainterPath()
        path.moveTo(0, 0)
        for i in range(self.num_lanes):
            cx = lane_w * i + lane_w/2
            path.lineTo(cx - well_width/2, 0)
            path.lineTo(cx - well_width/2, well_height)
            path.lineTo(cx + well_width/2, well_height)
            path.lineTo(cx + well_width/2, 0)
        path.lineTo(self.gel_width, 0)
        
        self.wells_path = QGraphicsPathItem(path, self)
        self.wells_path.setPen(QPen(Qt.GlobalColor.darkGray, 1))
        self.wells_path.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        for i in range(self.num_lanes):
            lane = GelLaneItem(self, index=i, width=lane_w, height=self.gel_height)
            lane.setPos(lane_w * i + lane_w/2, 0)
            lane.add_band(dist=40 + (i%2)*20)
            if i % 2 == 0:
                lane.add_band(dist=100 + i*10)
            self.lane_items.append(lane)

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, self.gel_width, self.gel_height)
        
    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(Qt.GlobalColor.blue, 2, Qt.PenStyle.DotLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect().adjusted(-2, -2, 2, 2))

    def to_dict(self):
        lanes_data = []
        for lane in self.lane_items:
            bands_data = []
            for band, _ in lane.bands:
                bands_data.append({
                    "y": band.y(),
                    "color": band.get_color().name(QColor.NameFormat.HexArgb)
                })
            lanes_data.append(bands_data)
            
        return {
            "type": "GelElectrophoresisItem",
            "pos": (self.pos().x(), self.pos().y()),
            "lanes": self.num_lanes,
            "width": self.gel_width,
            "height": self.gel_height,
            "show_labels": self.show_labels,
            "lanes_data": lanes_data
        }

    def load_dict(self, data):
        self.setPos(data.get("pos", (0,0))[0], data.get("pos", (0,0))[1])
        self.show_labels = data.get("show_labels", True)
        lanes_data = data.get("lanes_data", [])
        for i, lane_data in enumerate(lanes_data):
            if i < len(self.lane_items):
                lane = self.lane_items[i]
                for band, lbl in lane.bands:
                    self.scene().removeItem(band)
                    self.scene().removeItem(lbl)
                lane.bands.clear()
                
                for bdata in lane_data:
                    band = GelBandItem(lane, width=lane.lane_width*0.8)
                    band.setPos(0, bdata.get("y", lane._well_y + 40))
                    if "color" in bdata:
                        band.set_color(QColor(bdata["color"]))
                    
                    lbl = QGraphicsTextItem("", lane)
                    lbl.setFont(QFont("Arial", 8))
                    lane.bands.append((band, lbl))
                lane.update_labels()

    def to_json(self) -> dict:
        return self.to_dict()

    @classmethod
    def from_json(cls, data: dict) -> "PlateItem":
        item = GelElectrophoresisItem(lanes=data.get("lanes", 5), width=data.get("width", 250.0), height=data.get("height", 300.0))
        item.load_dict(data)
        return item

class PlateItem:
    @staticmethod
    def from_json(data: dict):
        t = data.get("type")
        if t == "TLCPlateItem":
            return TLCPlateItem.from_json(data)
        elif t == "GelElectrophoresisItem":
            return GelElectrophoresisItem.from_json(data)
        return None
