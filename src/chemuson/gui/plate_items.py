"""
Placas TLC y geles de electroforesis para Chemuson.

Manchas y bandas son items de nivel de escena (no hijos de la placa)
para que Qt entregue eventos del mouse correctamente.
"""
from __future__ import annotations

import math
from typing import Optional

from PyQt6.QtWidgets import (
    QGraphicsObject,
    QGraphicsItem,
    QGraphicsRectItem,
    QGraphicsEllipseItem,
    QGraphicsLineItem,
    QGraphicsTextItem,
    QGraphicsPathItem,
    QGraphicsSceneMouseEvent,
    QMenu,
)
from PyQt6.QtGui import (
    QColor,
    QPen,
    QBrush,
    QFont,
    QPainter,
    QPainterPath,
    QRadialGradient,
    QAction,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, pyqtSignal


# =============================================================================
# Gel scale conversion helpers (single source of truth)
# =============================================================================

def gel_normalized_migration(y: float, run_top_y: float, run_height: float) -> float:
    """Return t in [0,1]: 0=top (wells), 1=bottom (farthest migration)."""
    if run_height <= 0:
        return 0.0
    return max(0.0, min(1.0, (y - run_top_y) / run_height))


def gel_value_at_position(y: float, run_top_y: float, run_height: float,
                           scale_unit: str, mass_min_kda: float, mass_max_kda: float) -> float:
    """Convert a pixel Y position to the displayed value depending on scale mode.

    - Distance:       t  (0.00 top, 1.00 bottom)
    - log(Mass/kDa):  log10(mass_max) - t * (log10(mass_max) - log10(mass_min))
    - Mass(kDa):      10 ** (log10(mass_max) - t * (log10(mass_max) - log10(mass_min)))
    """
    t = gel_normalized_migration(y, run_top_y, run_height)
    log_min = math.log10(mass_min_kda) if mass_min_kda > 0 else 0.0
    log_max = math.log10(mass_max_kda) if mass_max_kda > 0 else 0.0

    if scale_unit == "Mass(kDa)":
        log_val = log_max - t * (log_max - log_min)
        return 10.0 ** log_val
    elif scale_unit == "log(Mass/kDa)":
        return log_max - t * (log_max - log_min)
    else:
        return t


def gel_position_for_value(value: float, run_top_y: float, run_height: float,
                            scale_unit: str, mass_min_kda: float, mass_max_kda: float) -> float:
    """Convert a displayed value to a pixel Y position (inverse of gel_value_at_position)."""
    log_min = math.log10(mass_min_kda) if mass_min_kda > 0 else 0.0
    log_max = math.log10(mass_max_kda) if mass_max_kda > 0 else 0.0

    if scale_unit == "Mass(kDa)":
        if value <= 0:
            return run_top_y + run_height
        log_val = math.log10(value)
        t = (log_max - log_val) / (log_max - log_min) if (log_max - log_min) != 0 else 0.0
    elif scale_unit == "log(Mass/kDa)":
        t = (log_max - value) / (log_max - log_min) if (log_max - log_min) != 0 else 0.0
    else:
        t = value

    t = max(0.0, min(1.0, t))
    return run_top_y + t * run_height


def gel_scale_ticks(scale_unit: str, mass_min_kda: float, mass_max_kda: float,
                    n_preferred: int = 7) -> list[float]:
    """Return nice tick values for the Y axis depending on scale mode.

    - Distance:       [0.00, 0.20, 0.40, 0.60, 0.80, 1.00]
    - log(Mass/kDa):  integer and half-integer log values within range
    - Mass(kDa):      decade + sub-decade values (1, 2, 5 pattern) within range
    """
    if scale_unit == "Distance":
        n = n_preferred - 1
        return [round(i / n, 2) for i in range(n + 1)]
    elif scale_unit == "log(Mass/kDa)":
        log_min = math.log10(mass_min_kda) if mass_min_kda > 0 else 0.0
        log_max = math.log10(mass_max_kda) if mass_max_kda > 0 else 0.0
        ticks = []
        step = 0.5
        v = math.ceil(log_max / step) * step
        while v >= log_min - 0.001:
            ticks.append(round(v, 2))
            v -= step
        return ticks
    else:
        log_min = math.log10(mass_min_kda) if mass_min_kda > 0 else 0.0
        log_max = math.log10(mass_max_kda) if mass_max_kda > 0 else 0.0
        sub_decade_multipliers = [1.0, 2.0, 5.0]
        ticks = []
        lo_exp = int(math.floor(log_min))
        hi_exp = int(math.ceil(log_max))
        for exp in range(hi_exp, lo_exp - 1, -1):
            for m in sub_decade_multipliers:
                val = m * (10.0 ** exp)
                if mass_min_kda - 0.001 <= val <= mass_max_kda + 0.001:
                    ticks.append(val)
        return ticks


# =============================================================================
# TLC Spot Item (scene-level)
# =============================================================================

class TLCSpotItem(QGraphicsObject):
    """Mancha individual en un carril TLC."""

    moved = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)

        self.spot_width = 18.0
        self.spot_height = 10.0
        self._color = QColor(80, 130, 190, 160)
        self.lane_ref: Optional[TLCLaneItem] = None
        self._show_rf_label = False

        self._setup_visual()

    def _setup_visual(self):
        self._gradient = QRadialGradient(0, 0, max(self.spot_width, self.spot_height) / 2)
        self._gradient.setColorAt(0, QColor(self._color.red(), self._color.green(), self._color.blue(), self._color.alpha()))
        self._gradient.setColorAt(0.6, QColor(self._color.red(), self._color.green(), self._color.blue(), int(self._color.alpha() * 0.6)))
        self._gradient.setColorAt(1.0, QColor(self._color.red(), self._color.green(), self._color.blue(), 0))

        self.spot_item = QGraphicsEllipseItem(
            -self.spot_width / 2, -self.spot_height / 2,
            self.spot_width, self.spot_height, self,
        )
        self.spot_item.setPen(QPen(Qt.PenStyle.NoPen))
        self.spot_item.setBrush(QBrush(self._gradient))

    def set_color(self, color: QColor):
        self._color = QColor(color)
        self._update_gradient()

    def get_color(self) -> QColor:
        return QColor(self._color)

    def set_opacity(self, opacity: float):
        alpha = max(0, min(255, int(opacity * 255)))
        self._color.setAlpha(alpha)
        self._update_gradient()

    def get_opacity(self) -> float:
        return self._color.alpha() / 255.0

    def _update_gradient(self):
        r, g, b = self._color.red(), self._color.green(), self._color.blue()
        a = self._color.alpha()
        self._gradient.setColorAt(0, QColor(r, g, b, a))
        self._gradient.setColorAt(0.6, QColor(r, g, b, int(a * 0.6)))
        self._gradient.setColorAt(1.0, QColor(r, g, b, 0))
        self.spot_item.setBrush(QBrush(self._gradient))

    def set_size(self, width: float, height: float):
        self.spot_width = max(4.0, width)
        self.spot_height = max(2.0, height)
        self.spot_item.setRect(
            -self.spot_width / 2, -self.spot_height / 2,
            self.spot_width, self.spot_height,
        )
        self._gradient.setRadius(max(self.spot_width, self.spot_height) / 2)
        self._update_gradient()

    def boundingRect(self) -> QRectF:
        return self.spot_item.boundingRect().adjusted(-1, -1, 1, 1)

    def shape(self) -> QPainterPath:
        path = QPainterPath()
        path.addEllipse(
            -self.spot_width / 2, -self.spot_height / 2,
            self.spot_width, self.spot_height,
        )
        return path

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(59, 130, 246), 1.5, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawEllipse(
                int(-self.spot_width / 2 - 2), int(-self.spot_height / 2 - 2),
                int(self.spot_width + 4), int(self.spot_height + 4),
            )

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange and self.lane_ref is not None:
            new_pos = QPointF(value)
            new_pos.setX(self.lane_ref.scenePos().x())
            y = new_pos.y()
            min_y = self.lane_ref.scenePos().y() + self.lane_ref.solvent_front_y()
            max_y = self.lane_ref.scenePos().y() + self.lane_ref.baseline_y()
            y = max(min_y, min(max_y, y))
            new_pos.setY(y)
            return new_pos
        elif change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.moved.emit()
            if self.lane_ref is not None:
                self.lane_ref.update_rf_labels()
        return super().itemChange(change, value)

    def contextMenuEvent(self, event: QGraphicsSceneMouseEvent):
        menu = QMenu()
        color_action = QAction("Cambiar color...", menu)
        color_action.triggered.connect(self._change_color)
        menu.addAction(color_action)

        opacity_menu = menu.addMenu("Opacidad")
        for pct, alpha in [("20%", 51), ("40%", 102), ("60%", 153), ("80%", 204), ("100%", 255)]:
            act = QAction(pct, menu)
            act.triggered.connect(lambda checked, a=alpha: self.set_opacity(a / 255.0))
            opacity_menu.addAction(act)

        size_menu = menu.addMenu("Tamano")
        for label, w, h in [("Pequena", 12, 6), ("Normal", 18, 10), ("Grande", 26, 14)]:
            act = QAction(label, menu)
            act.triggered.connect(lambda checked, ww=w, hh=h: self.set_size(ww, hh))
            size_menu.addAction(act)

        show_rf_action = QAction("Mostrar Rf" if not self._show_rf_label else "Ocultar Rf", menu)
        show_rf_action.triggered.connect(self._toggle_rf_label)
        menu.addAction(show_rf_action)

        set_rf_action = QAction("Set Rf...", menu)
        set_rf_action.triggered.connect(self._set_rf)
        menu.addAction(set_rf_action)

        delete_action = QAction("Eliminar mancha", menu)
        delete_action.triggered.connect(self._delete_self)
        menu.addAction(delete_action)

        menu.exec(event.screenPos())

    def _toggle_rf_label(self):
        self._show_rf_label = not self._show_rf_label
        if self.lane_ref is not None:
            self.lane_ref.update_rf_labels()

    def _set_rf(self):
        from PyQt6.QtWidgets import QInputDialog
        if self.lane_ref is None:
            return
        plate = self.lane_ref.parentItem()
        if plate is None:
            return
        dist_total = self.lane_ref.baseline_y() - self.lane_ref.solvent_front_y()
        dist_spot = self.lane_ref.baseline_y() - (self.y() - self.lane_ref.scenePos().y())
        current_rf = dist_spot / dist_total if dist_total != 0 else 0.0
        current_rf = max(0.0, min(1.0, current_rf))
        rf_val, ok = QInputDialog.getDouble(
            None, "Set Rf", "Valor de Rf (0.0 - 1.0):",
            current_rf, 0.0, 1.0, 3,
        )
        if not ok:
            return
        rf_val = max(0.0, min(1.0, rf_val))
        target_y = self.lane_ref.scenePos().y() + self.lane_ref.baseline_y() - rf_val * dist_total
        self.setY(target_y)
        self._show_rf_label = True
        if self.lane_ref is not None:
            self.lane_ref.update_rf_labels()

    def _change_color(self):
        from PyQt6.QtWidgets import QColorDialog
        picked = QColorDialog.getColor(self._color, None, "Color de mancha")
        if picked.isValid():
            self.set_color(picked)

    def _delete_self(self):
        if self.lane_ref is not None:
            self.lane_ref.remove_spot(self)

    def to_dict(self) -> dict:
        return {
            "x": self.x(),
            "y": self.y(),
            "color": self._color.name(QColor.NameFormat.HexArgb),
            "width": self.spot_width,
            "height": self.spot_height,
        }

    def load_dict(self, data: dict):
        self.setPos(data.get("x", 0), data.get("y", 0))
        if "color" in data:
            self.set_color(QColor(data["color"]))
        if "width" in data and "height" in data:
            self.set_size(data["width"], data["height"])


# =============================================================================
# TLC Lane Item
# =============================================================================

class TLCLaneItem(QGraphicsObject):
    """Carril individual en una placa TLC."""

    def __init__(self, parent=None, index=0, width=40.0, height=200.0):
        super().__init__(parent)
        self.index = index
        self.lane_width = width
        self.lane_height = height
        self._baseline_y = height * 0.85
        self._solvent_front_y = height * 0.15
        self.rf_labels: list[tuple[TLCSpotItem, QGraphicsTextItem]] = []

    def set_boundaries(self, baseline_y: float, solvent_front_y: float):
        self._baseline_y = baseline_y
        self._solvent_front_y = solvent_front_y
        self.update_rf_labels()

    def baseline_y(self) -> float:
        return self._baseline_y

    def solvent_front_y(self) -> float:
        return self._solvent_front_y

    def add_spot(self, rf: float = 0.5, scene=None) -> TLCSpotItem:
        spot = TLCSpotItem()
        spot.lane_ref = self
        cx = self.scenePos().x()
        y = self.scenePos().y() + self._baseline_y - rf * (self._baseline_y - self._solvent_front_y)
        spot.setPos(cx, y)
        if scene is not None:
            scene.addItem(spot)

        label = QGraphicsTextItem("", self)
        font = QFont("Arial", 8)
        label.setFont(font)
        label.setDefaultTextColor(QColor(30, 41, 59))
        self.rf_labels.append((spot, label))
        self.update_rf_labels()
        return spot

    def remove_spot(self, spot: TLCSpotItem):
        for i, (s, lbl) in enumerate(self.rf_labels):
            if s is spot:
                if self.scene():
                    self.scene().removeItem(s)
                    self.scene().removeItem(lbl)
                del self.rf_labels[i]
                self.update_rf_labels()
                break

    def update_rf_labels(self):
        plate = self.parentItem()
        show_labels = getattr(plate, "show_labels", True)
        dist_total = self._baseline_y - self._solvent_front_y

        for spot, label in self.rf_labels:
            if not show_labels or not getattr(spot, "_show_rf_label", False):
                label.setVisible(False)
                continue
            label.setVisible(True)
            dist_spot = self._baseline_y - (spot.y() - self.scenePos().y())
            rf = dist_spot / dist_total if dist_total != 0 else 0.0
            rf = max(0.0, min(1.0, rf))
            label.setPlainText(f"{rf:.2f}")
            lbl_w = label.boundingRect().width()
            label.setPos(-lbl_w / 2, spot.y() - self.scenePos().y() - label.boundingRect().height() - 2)

    def boundingRect(self) -> QRectF:
        return QRectF(-self.lane_width / 2, 0, self.lane_width, self.lane_height)

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(100, 100, 255, 100), 1, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect())

    def to_dict(self) -> dict:
        spots = []
        for spot, _ in self.rf_labels:
            spots.append(spot.to_dict())
        return {"spots": spots}

    def load_dict(self, data: dict, scene=None):
        for spot, lbl in self.rf_labels:
            if self.scene():
                self.scene().removeItem(spot)
                self.scene().removeItem(lbl)
        self.rf_labels.clear()

        for sdata in data.get("spots", []):
            spot = TLCSpotItem()
            spot.lane_ref = self
            spot.load_dict(sdata)
            if scene is not None:
                scene.addItem(spot)
            lbl = QGraphicsTextItem("", self)
            lbl.setFont(QFont("Arial", 8))
            lbl.setDefaultTextColor(QColor(30, 41, 59))
            self.rf_labels.append((spot, lbl))
        self.update_rf_labels()


# =============================================================================
# TLC Plate Item
# =============================================================================

class TLCPlateItem(QGraphicsObject):
    """Placa de cromatografia en capa fina (TLC)."""

    def __init__(self, lanes: int = 3, width: float = 180.0, height: float = 220.0, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)

        self.num_lanes = lanes
        self.plate_width = max(width, lanes * 36.0)
        self.plate_height = height
        self.baseline_y = self.plate_height * 0.85
        self.solvent_front_y = self.plate_height * 0.15
        self.show_labels = True
        self.lane_items: list[TLCLaneItem] = []
        self._last_plate_pos = QPointF(0, 0)

        self._setup_plate()

    def _setup_plate(self):
        self.outline = QGraphicsRectItem(0, 0, self.plate_width, self.plate_height, self)
        self.outline.setPen(QPen(QColor(30, 41, 59), 1.5))
        self.outline.setBrush(QBrush(QColor(255, 255, 255, 240)))

        pen_dash = QPen(QColor(100, 116, 139), 1, Qt.PenStyle.DashLine)
        self.baseline_line = QGraphicsLineItem(
            8, self.baseline_y, self.plate_width - 8, self.baseline_y, self,
        )
        self.baseline_line.setPen(pen_dash)

        self.solvent_line = QGraphicsLineItem(
            8, self.solvent_front_y, self.plate_width - 8, self.solvent_front_y, self,
        )
        self.solvent_line.setPen(pen_dash)

        lane_w = self.plate_width / self.num_lanes
        for i in range(self.num_lanes):
            lane = TLCLaneItem(self, index=i, width=lane_w, height=self.plate_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            lane.set_boundaries(self.baseline_y, self.solvent_front_y)
            self.lane_items.append(lane)

    def add_spots_to_lanes(self, scene=None, rf_values=None):
        for i, lane in enumerate(self.lane_items):
            rf = rf_values[i] if rf_values and i < len(rf_values) else 0.0
            lane.add_spot(rf=rf, scene=scene)

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, self.plate_width, self.plate_height)

    def shape(self) -> QPainterPath:
        path = QPainterPath()
        path.addRect(self.boundingRect())
        return path

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(59, 130, 246), 2, Qt.PenStyle.DotLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect().adjusted(-3, -3, 3, 3))

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            new_pos = QPointF(value)
            delta = new_pos - self._last_plate_pos
            for lane in self.lane_items:
                for spot, _ in lane.rf_labels:
                    spot.setPos(spot.pos() + delta)
            self._last_plate_pos = new_pos
        return super().itemChange(change, value)

    def set_num_lanes(self, new_count: int, scene=None):
        if new_count == self.num_lanes:
            return
        old_data = [lane.to_dict() for lane in self.lane_items]
        for lane in self.lane_items:
            for spot, lbl in lane.rf_labels:
                if self.scene():
                    self.scene().removeItem(spot)
                    self.scene().removeItem(lbl)
            if self.scene():
                self.scene().removeItem(lane)
        self.lane_items.clear()
        self.num_lanes = new_count
        self.plate_width = max(self.plate_width, new_count * 36.0)
        self.outline.setRect(0, 0, self.plate_width, self.plate_height)
        self.baseline_line.setLine(8, self.baseline_y, self.plate_width - 8, self.baseline_y)
        self.solvent_line.setLine(8, self.solvent_front_y, self.plate_width - 8, self.solvent_front_y)

        lane_w = self.plate_width / self.num_lanes
        for i in range(self.num_lanes):
            lane = TLCLaneItem(self, index=i, width=lane_w, height=self.plate_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            lane.set_boundaries(self.baseline_y, self.solvent_front_y)
            if i < len(old_data):
                lane.load_dict(old_data[i], scene=scene)
            else:
                lane.add_spot(rf=0.0, scene=scene)
            self.lane_items.append(lane)
        self._last_plate_pos = self.pos()

    def contextMenuEvent(self, event: QGraphicsSceneMouseEvent):
        menu = QMenu()

        lanes_menu = menu.addMenu("Numero de carriles")
        for n in range(1, 13):
            act = QAction(str(n), menu)
            act.triggered.connect(lambda checked, nn=n: self.set_num_lanes(nn, scene=self.scene()))
            lanes_menu.addAction(act)

        label_action = QAction("Mostrar etiquetas" if not self.show_labels else "Ocultar etiquetas", menu)
        label_action.triggered.connect(self._toggle_labels)
        menu.addAction(label_action)

        add_lane_action = QAction("Anadir carril", menu)
        add_lane_action.triggered.connect(lambda: self.set_num_lanes(self.num_lanes + 1, scene=self.scene()))
        menu.addAction(add_lane_action)

        add_spot_menu = menu.addMenu("Add spot")
        for i, lane in enumerate(self.lane_items):
            act = QAction(f"Carril {i + 1}", menu)
            act.triggered.connect(lambda checked, ln=lane: self._add_spot_to_lane(ln))
            add_spot_menu.addAction(act)

        add_line_action = QAction("Add line", menu)
        add_line_action.triggered.connect(self._add_extra_line)
        menu.addAction(add_line_action)

        menu.exec(event.screenPos())

    def _add_spot_to_lane(self, lane):
        if self.scene() is None:
            return
        lane.add_spot(rf=0.5, scene=self.scene())

    def _add_extra_line(self):
        from PyQt6.QtWidgets import QInputDialog
        options = ["Baseline", "Solvent front", "Custom line"]
        line_type, ok = QInputDialog.getItem(
            None, "Add line", "Tipo de linea:", options, 2, False,
        )
        if not ok:
            return
        if line_type == "Baseline":
            y_pos, ok = QInputDialog.getDouble(
                None, "Baseline Y", "Posicion Y (% del alto):",
                85.0, 10.0, 95.0, 1,
            )
            if ok:
                self.baseline_y = self.plate_height * y_pos / 100.0
                self.baseline_line.setLine(8, self.baseline_y, self.plate_width - 8, self.baseline_y)
                for lane in self.lane_items:
                    lane.set_boundaries(self.baseline_y, self.solvent_front_y)
        elif line_type == "Solvent front":
            y_pos, ok = QInputDialog.getDouble(
                None, "Solvent front Y", "Posicion Y (% del alto):",
                15.0, 5.0, 90.0, 1,
            )
            if ok:
                self.solvent_front_y = self.plate_height * y_pos / 100.0
                self.solvent_line.setLine(8, self.solvent_front_y, self.plate_width - 8, self.solvent_front_y)
                for lane in self.lane_items:
                    lane.set_boundaries(self.baseline_y, self.solvent_front_y)
        else:
            y_pos, ok = QInputDialog.getDouble(
                None, "Custom line Y", "Posicion Y (% del alto):",
                50.0, 0.0, 100.0, 1,
            )
            if ok:
                pen_dash = QPen(QColor(100, 116, 139), 1, Qt.PenStyle.DotLine)
                extra_line = QGraphicsLineItem(
                    8, y_pos, self.plate_width - 8, y_pos, self,
                )
                extra_line.setPen(pen_dash)

    def _toggle_labels(self):
        self.show_labels = not self.show_labels
        for lane in self.lane_items:
            lane.update_rf_labels()

    def to_dict(self) -> dict:
        lanes_data = []
        for lane in self.lane_items:
            lanes_data.append(lane.to_dict())
        return {
            "type": "TLCPlateItem",
            "pos": (self.pos().x(), self.pos().y()),
            "lanes": self.num_lanes,
            "width": self.plate_width,
            "height": self.plate_height,
            "show_labels": self.show_labels,
            "lanes_data": lanes_data,
        }

    def load_dict(self, data: dict, scene=None):
        self.setPos(data.get("pos", (0, 0))[0], data.get("pos", (0, 0))[1])
        self.show_labels = data.get("show_labels", True)
        self.num_lanes = data.get("lanes", 3)
        self.plate_width = data.get("width", 180.0)
        self.plate_height = data.get("height", 220.0)
        self.baseline_y = self.plate_height * 0.85
        self.solvent_front_y = self.plate_height * 0.15

        for lane in self.lane_items:
            for spot, lbl in lane.rf_labels:
                if self.scene():
                    self.scene().removeItem(spot)
                    self.scene().removeItem(lbl)
            if self.scene():
                self.scene().removeItem(lane)
        self.lane_items.clear()

        self.outline.setRect(0, 0, self.plate_width, self.plate_height)
        self.baseline_line.setLine(8, self.baseline_y, self.plate_width - 8, self.baseline_y)
        self.solvent_line.setLine(8, self.solvent_front_y, self.plate_width - 8, self.solvent_front_y)

        lanes_data = data.get("lanes_data", [])
        lane_w = self.plate_width / self.num_lanes
        for i in range(self.num_lanes):
            lane = TLCLaneItem(self, index=i, width=lane_w, height=self.plate_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            lane.set_boundaries(self.baseline_y, self.solvent_front_y)
            if i < len(lanes_data):
                lane.load_dict(lanes_data[i], scene=scene)
            else:
                lane.add_spot(rf=0.0, scene=scene)
            self.lane_items.append(lane)

    def to_json(self) -> dict:
        return self.to_dict()

    @classmethod
    def from_json(cls, data: dict) -> "TLCPlateItem":
        item = cls(
            lanes=data.get("lanes", 3),
            width=data.get("width", 180.0),
            height=data.get("height", 220.0),
        )
        item.load_dict(data)
        return item


# =============================================================================
# Gel Band Item (scene-level)
# =============================================================================

class GelBandItem(QGraphicsObject):
    """Banda individual en un carril de gel de electroforesis."""

    moved = pyqtSignal()

    def __init__(self, parent=None, width=32.0):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)

        self.band_width = width
        self.band_height = 7.0
        self._color = QColor(40, 40, 50, 180)
        self.lane_ref: Optional[GelLaneItem] = None
        self._show_label = False

        self._setup_visual()

    def _setup_visual(self):
        self.band_rect = QGraphicsRectItem(
            -self.band_width / 2, -self.band_height / 2,
            self.band_width, self.band_height, self,
        )
        self.band_rect.setPen(QPen(Qt.PenStyle.NoPen))
        self.band_rect.setBrush(QBrush(self._color))

    def set_color(self, color: QColor):
        self._color = QColor(color)
        self._update_brush()

    def get_color(self) -> QColor:
        return QColor(self._color)

    def set_opacity(self, opacity: float):
        alpha = max(0, min(255, int(opacity * 255)))
        self._color.setAlpha(alpha)
        self._update_brush()

    def get_opacity(self) -> float:
        return self._color.alpha() / 255.0

    def _update_brush(self):
        self.band_rect.setBrush(QBrush(self._color))

    def set_size(self, width: float, height: float):
        self.band_width = max(4.0, width)
        self.band_height = max(2.0, height)
        self.band_rect.setRect(
            -self.band_width / 2, -self.band_height / 2,
            self.band_width, self.band_height,
        )

    def boundingRect(self) -> QRectF:
        return self.band_rect.boundingRect().adjusted(-1, -1, 1, 1)

    def shape(self) -> QPainterPath:
        path = QPainterPath()
        path.addRect(
            -self.band_width / 2, -self.band_height / 2,
            self.band_width, self.band_height,
        )
        return path

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(59, 130, 246), 1.5, Qt.PenStyle.DashLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(
                int(-self.band_width / 2 - 2), int(-self.band_height / 2 - 2),
                int(self.band_width + 4), int(self.band_height + 4),
            )

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionChange and self.lane_ref is not None:
            new_pos = QPointF(value)
            new_pos.setX(self.lane_ref.scenePos().x())
            y = new_pos.y()
            min_y = self.lane_ref.scenePos().y() + self.lane_ref.well_y()
            max_y = self.lane_ref.scenePos().y() + self.lane_ref.lane_height - 5
            y = max(min_y, min(max_y, y))
            new_pos.setY(y)
            return new_pos
        elif change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.moved.emit()
            if self.lane_ref is not None:
                self.lane_ref.update_labels()
        return super().itemChange(change, value)

    def contextMenuEvent(self, event: QGraphicsSceneMouseEvent):
        menu = QMenu()

        color_action = QAction("Cambiar color...", menu)
        color_action.triggered.connect(self._change_color)
        menu.addAction(color_action)

        opacity_menu = menu.addMenu("Opacidad")
        for pct, alpha in [("20%", 51), ("40%", 102), ("60%", 153), ("80%", 204), ("100%", 255)]:
            act = QAction(pct, menu)
            act.triggered.connect(lambda checked, a=alpha: self.set_opacity(a / 255.0))
            opacity_menu.addAction(act)

        size_menu = menu.addMenu("Ancho de banda")
        for label, w in [("Estrecha", 20), ("Normal", 32), ("Ancha", 44)]:
            act = QAction(label, menu)
            act.triggered.connect(lambda checked, ww=w: self.set_size(ww, self.band_height))
            size_menu.addAction(act)

        show_label_action = QAction("Mostrar etiqueta" if not self._show_label else "Ocultar etiqueta", menu)
        show_label_action.triggered.connect(self._toggle_label)
        menu.addAction(show_label_action)

        set_dist_action = QAction("Set distancia...", menu)
        set_dist_action.triggered.connect(self._set_distance)
        menu.addAction(set_dist_action)

        delete_action = QAction("Eliminar banda", menu)
        delete_action.triggered.connect(self._delete_self)
        menu.addAction(delete_action)

        menu.exec(event.screenPos())

    def _toggle_label(self):
        self._show_label = not self._show_label
        if self.lane_ref is not None:
            self.lane_ref.update_labels()

    def _set_distance(self):
        from PyQt6.QtWidgets import QInputDialog
        if self.lane_ref is None:
            return
        gel = self.lane_ref.parentItem()
        if gel is None:
            return
        dist = self.lane_ref._well_y - (self.y() - self.lane_ref.scenePos().y())
        current_cm = dist * 0.05
        dist_val, ok = QInputDialog.getDouble(
            None, "Set distancia", "Distancia migrada (cm):",
            current_cm, 0.0, 50.0, 2,
        )
        if not ok:
            return
        target_dist = dist_val / 0.05
        target_y = self.lane_ref.scenePos().y() + self.lane_ref._well_y - target_dist
        min_y = self.lane_ref.scenePos().y() + self.lane_ref.lane_height - 5
        max_y = self.lane_ref.scenePos().y() + self.lane_ref.well_y()
        target_y = max(min_y, min(max_y, target_y))
        self.setY(target_y)
        self._show_label = True
        if self.lane_ref is not None:
            self.lane_ref.update_labels()

    def _change_color(self):
        from PyQt6.QtWidgets import QColorDialog
        picked = QColorDialog.getColor(self._color, None, "Color de banda")
        if picked.isValid():
            self.set_color(picked)

    def _delete_self(self):
        if self.lane_ref is not None:
            self.lane_ref.remove_band(self)

    def to_dict(self) -> dict:
        return {
            "x": self.x(),
            "y": self.y(),
            "color": self._color.name(QColor.NameFormat.HexArgb),
            "width": self.band_width,
            "height": self.band_height,
        }

    def load_dict(self, data: dict):
        self.setPos(data.get("x", 0), data.get("y", 0))
        if "color" in data:
            self.set_color(QColor(data["color"]))
        if "width" in data and "height" in data:
            self.set_size(data["width"], data["height"])


# =============================================================================
# Gel Lane Item
# =============================================================================

class GelLaneItem(QGraphicsObject):
    """Carril individual en un gel de electroforesis.

    Los pozos estan en la parte superior; las bandas migran hacia abajo.
    """

    def __init__(self, parent=None, index=0, width=40.0, height=250.0):
        super().__init__(parent)
        self.index = index
        self.lane_width = width
        self.lane_height = height
        self._well_y = 22.0
        self.bands: list[tuple[GelBandItem, QGraphicsTextItem]] = []

    def well_y(self) -> float:
        return self._well_y

    def add_band(self, dist: float = 0.0, scene=None) -> GelBandItem:
        band = GelBandItem(width=self.lane_width * 0.75)
        band.lane_ref = self
        cx = self.scenePos().x()
        y = self.scenePos().y() + self._well_y + dist
        band.setPos(cx, y)
        if scene is not None:
            scene.addItem(band)

        label = QGraphicsTextItem("", self)
        font = QFont("Arial", 8)
        label.setFont(font)
        label.setDefaultTextColor(QColor(30, 41, 59))
        self.bands.append((band, label))
        self.update_labels()
        return band

    def remove_band(self, band: GelBandItem):
        for i, (b, lbl) in enumerate(self.bands):
            if b is band:
                if self.scene():
                    self.scene().removeItem(b)
                    self.scene().removeItem(lbl)
                del self.bands[i]
                self.update_labels()
                break

    def update_labels(self):
        gel = self.parentItem()
        show_labels = getattr(gel, "show_labels", True)
        scale_unit = getattr(gel, "scale_unit", "Distance")
        mass_min_kda = getattr(gel, "mass_min_kda", 10.0)
        mass_max_kda = getattr(gel, "mass_max_kda", 250.0)
        run_top_y = self._well_y
        run_height = self.lane_height - self._well_y

        for band, label in self.bands:
            if not show_labels or not getattr(band, "_show_label", False):
                label.setVisible(False)
                continue
            label.setVisible(True)
            band_y_local = (band.y() - self.scenePos().y())
            val = gel_value_at_position(band_y_local, run_top_y, run_height,
                                         scale_unit, mass_min_kda, mass_max_kda)
            if scale_unit == "Distance":
                label.setPlainText(f"{val:.2f}")
            elif scale_unit == "log(Mass/kDa)":
                label.setPlainText(f"{val:.2f}")
            else:
                label.setPlainText(f"{val:.0f}")
            lbl_w = label.boundingRect().width()
            label.setPos(-lbl_w / 2, band.y() - self.scenePos().y() - label.boundingRect().height() - 2)

    def boundingRect(self) -> QRectF:
        return QRectF(-self.lane_width / 2, 0, self.lane_width, self.lane_height)

    def paint(self, painter, option, widget=None):
        pass

    def to_dict(self) -> dict:
        bands_data = []
        for band, _ in self.bands:
            bands_data.append(band.to_dict())
        return {"bands": bands_data}

    def load_dict(self, data: dict, scene=None):
        for band, lbl in self.bands:
            if self.scene():
                self.scene().removeItem(band)
                self.scene().removeItem(lbl)
        self.bands.clear()

        for bdata in data.get("bands", []):
            band = GelBandItem(width=self.lane_width * 0.75)
            band.lane_ref = self
            band.load_dict(bdata)
            if scene is not None:
                scene.addItem(band)
            lbl = QGraphicsTextItem("", self)
            lbl.setFont(QFont("Arial", 8))
            lbl.setDefaultTextColor(QColor(30, 41, 59))
            self.bands.append((band, lbl))
        self.update_labels()


# =============================================================================
# Gel Electrophoresis Item
# =============================================================================

class GelElectrophoresisItem(QGraphicsObject):
    """Gel de electroforesis con pozos y bandas.

    Los pozos estan en la parte superior; las bandas migran hacia abajo.
    """

    def __init__(self, lanes: int = 5, width: float = 280.0, height: float = 320.0, parent=None):
        super().__init__(parent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)

        self.num_lanes = lanes
        self.gel_width = max(width, lanes * 44.0)
        self.gel_height = height
        self.show_labels = True
        self.show_scale = True
        self.scale_unit = "Distance"
        self.scale_min = 0.0
        self.scale_max = 1.0
        self.mass_min_kda = 10.0
        self.mass_max_kda = 250.0
        self.lane_items: list[GelLaneItem] = []
        self._last_gel_pos = QPointF(0, 0)
        self._scale_lines: list[QGraphicsLineItem] = []
        self._scale_labels: list[QGraphicsTextItem] = []
        self._scale_title: Optional[QGraphicsTextItem] = None

        self._setup_gel()

    def _setup_gel(self):
        self.outline = QGraphicsRectItem(0, 0, self.gel_width, self.gel_height, self)
        self.outline.setPen(QPen(QColor(71, 85, 105), 1.5))
        self.outline.setBrush(QBrush(QColor(255, 255, 255, 255)))

        lane_w = self.gel_width / self.num_lanes
        well_height = 14.0
        well_width = lane_w * 0.65

        path = QPainterPath()
        path.moveTo(0, 0)
        for i in range(self.num_lanes):
            cx = lane_w * i + lane_w / 2
            path.lineTo(cx - well_width / 2, 0)
            path.lineTo(cx - well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, 0)
        path.lineTo(self.gel_width, 0)
        path.lineTo(self.gel_width, self.gel_height)
        path.lineTo(0, self.gel_height)
        path.closeSubpath()

        self.wells_path = QGraphicsPathItem(path, self)
        self.wells_path.setPen(QPen(QColor(71, 85, 105), 1.5))
        self.wells_path.setBrush(QBrush(QColor(255, 255, 255, 255)))

        for i in range(self.num_lanes):
            lane = GelLaneItem(self, index=i, width=lane_w, height=self.gel_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            self.lane_items.append(lane)

        self._build_scale()

    def _build_scale(self):
        for line in self._scale_lines:
            if line.scene():
                line.scene().removeItem(line)
        for label in self._scale_labels:
            if label.scene():
                label.scene().removeItem(label)
        if self._scale_title is not None and self._scale_title.scene():
            self._scale_title.scene().removeItem(self._scale_title)
        self._scale_lines.clear()
        self._scale_labels.clear()
        self._scale_title = None

        if not self.show_scale:
            return

        pen = QPen(QColor(100, 116, 139), 1, Qt.PenStyle.DashLine)
        font = QFont("Arial", 8)
        gel_left = 0
        run_top_y = 22.0
        run_bottom_y = self.gel_height
        run_height = run_bottom_y - run_top_y

        tick_values = gel_scale_ticks(self.scale_unit, self.mass_min_kda, self.mass_max_kda)

        for val in tick_values:
            y_local = gel_position_for_value(val, run_top_y, run_height,
                                              self.scale_unit, self.mass_min_kda, self.mass_max_kda)
            line = QGraphicsLineItem(gel_left - 14, y_local, gel_left - 2, y_local)
            line.setPen(pen)
            line.setParentItem(self)
            self._scale_lines.append(line)

            if self.scale_unit == "Distance":
                text = f"{val:.2f}"
            elif self.scale_unit == "log(Mass/kDa)":
                text = f"{val:.2f}"
            else:
                text = f"{val:.0f}"
            label = QGraphicsTextItem(text)
            label.setFont(font)
            label.setDefaultTextColor(QColor(30, 41, 59))
            lbl_w = label.boundingRect().width()
            lbl_h = label.boundingRect().height()
            label.setPos(gel_left - 16 - lbl_w, y_local - lbl_h / 2)
            label.setParentItem(self)
            self._scale_labels.append(label)

        title_font = QFont("Arial", 8)
        title_font.setBold(True)
        title_text = self.scale_unit
        title_label = QGraphicsTextItem(title_text)
        title_label.setFont(title_font)
        title_label.setDefaultTextColor(QColor(30, 41, 59))
        title_lbl_w = title_label.boundingRect().width()
        title_lbl_h = title_label.boundingRect().height()
        mid_y = run_top_y + run_height / 2
        max_label_w = max((t.boundingRect().width() for t in self._scale_labels), default=0)
        center_x = gel_left - 16 - max_label_w - title_lbl_h / 2 - 10
        center_y = mid_y
        title_label.setPos(center_x - title_lbl_w / 2, center_y - title_lbl_h / 2)
        title_label.setTransformOriginPoint(title_lbl_w / 2, title_lbl_h / 2)
        title_label.setRotation(-90)
        title_label.setParentItem(self)
        self._scale_title = title_label

    def _well_y(self):
        return 22.0

    def add_bands_to_lanes(self, scene=None, dist_values=None):
        for i, lane in enumerate(self.lane_items):
            dist = dist_values[i] if dist_values and i < len(dist_values) else 0.0
            lane.add_band(dist=dist, scene=scene)

    def boundingRect(self) -> QRectF:
        return QRectF(0, 0, self.gel_width, self.gel_height)

    def shape(self) -> QPainterPath:
        path = QPainterPath()
        path.addRect(self.boundingRect())
        return path

    def paint(self, painter, option, widget=None):
        if self.isSelected():
            painter.setPen(QPen(QColor(59, 130, 246), 2, Qt.PenStyle.DotLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(self.boundingRect().adjusted(-3, -3, 3, 3))

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            new_pos = QPointF(value)
            delta = new_pos - self._last_gel_pos
            for lane in self.lane_items:
                for band, _ in lane.bands:
                    band.setPos(band.pos() + delta)
            self._last_gel_pos = new_pos
        return super().itemChange(change, value)

    def set_num_lanes(self, new_count: int, scene=None):
        if new_count == self.num_lanes:
            return
        old_data = [lane.to_dict() for lane in self.lane_items]
        for lane in self.lane_items:
            for band, lbl in lane.bands:
                if self.scene():
                    self.scene().removeItem(band)
                    self.scene().removeItem(lbl)
            if self.scene():
                self.scene().removeItem(lane)
        self.lane_items.clear()
        self.num_lanes = new_count
        self.gel_width = max(self.gel_width, new_count * 44.0)

        lane_w = self.gel_width / self.num_lanes
        well_height = 14.0
        well_width = lane_w * 0.65

        path = QPainterPath()
        path.moveTo(0, 0)
        for i in range(self.num_lanes):
            cx = lane_w * i + lane_w / 2
            path.lineTo(cx - well_width / 2, 0)
            path.lineTo(cx - well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, 0)
        path.lineTo(self.gel_width, 0)
        path.lineTo(self.gel_width, self.gel_height)
        path.lineTo(0, self.gel_height)
        path.closeSubpath()

        self.outline.setRect(0, 0, self.gel_width, self.gel_height)
        self.wells_path.setPath(path)

        for i in range(self.num_lanes):
            lane = GelLaneItem(self, index=i, width=lane_w, height=self.gel_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            if i < len(old_data):
                lane.load_dict(old_data[i], scene=scene)
            else:
                lane.add_band(dist=0.0, scene=scene)
            self.lane_items.append(lane)
        self._build_scale()

    def contextMenuEvent(self, event: QGraphicsSceneMouseEvent):
        menu = QMenu()

        lanes_menu = menu.addMenu("Numero de carriles")
        for n in range(1, 21):
            act = QAction(str(n), menu)
            act.triggered.connect(lambda checked, nn=n: self.set_num_lanes(nn, scene=self.scene()))
            lanes_menu.addAction(act)

        label_action = QAction("Mostrar etiquetas" if not self.show_labels else "Ocultar etiquetas", menu)
        label_action.triggered.connect(self._toggle_labels)
        menu.addAction(label_action)

        scale_action = QAction("Mostrar escala" if not self.show_scale else "Ocultar escala", menu)
        scale_action.triggered.connect(self._toggle_scale)
        menu.addAction(scale_action)

        add_lane_action = QAction("Anadir carril", menu)
        add_lane_action.triggered.connect(lambda: self.set_num_lanes(self.num_lanes + 1, scene=self.scene()))
        menu.addAction(add_lane_action)

        add_band_menu = menu.addMenu("Add band")
        for i, lane in enumerate(self.lane_items):
            act = QAction(f"Carril {i + 1}", menu)
            act.triggered.connect(lambda checked, ln=lane: self._add_band_to_lane(ln))
            add_band_menu.addAction(act)

        scale_unit_menu = menu.addMenu("Unidad de escala")
        for unit in ["Distance", "Mass(kDa)", "log(Mass/kDa)"]:
            act = QAction(unit, menu)
            act.setCheckable(True)
            act.setChecked(self.scale_unit == unit)
            act.triggered.connect(lambda checked, u=unit: self._set_scale_unit(u))
            scale_unit_menu.addAction(act)

        range_action = QAction("Configurar rango de escala...", menu)
        range_action.triggered.connect(self._set_scale_range)
        menu.addAction(range_action)

        menu.exec(event.screenPos())

    def _add_band_to_lane(self, lane):
        if self.scene() is None:
            return
        lane.add_band(dist=0.0, scene=self.scene())

    def _set_scale_unit(self, unit):
        self.scale_unit = unit
        self._build_scale()
        for lane in self.lane_items:
            lane.update_labels()

    def _set_scale_range(self):
        from PyQt6.QtWidgets import QInputDialog
        if self.scale_unit == "Distance":
            mn, ok1 = QInputDialog.getDouble(
                None, "Escala minima", "Valor minimo (0.0 para Distance):",
                self.scale_min, 0.0, 1.0, 2,
            )
            if not ok1:
                return
            mx, ok2 = QInputDialog.getDouble(
                None, "Escala maxima", "Valor maximo (1.0 para Distance):",
                self.scale_max, mn, 1.0, 2,
            )
            if not ok2:
                return
            self.scale_min = mn
            self.scale_max = mx
        elif self.scale_unit == "Mass(kDa)":
            mn, ok1 = QInputDialog.getDouble(
                None, "Masa minima (kDa)", "Masa minima en kDa (banda mas pequena):",
                self.mass_min_kda, 0.1, 1000000.0, 2,
            )
            if not ok1:
                return
            mx, ok2 = QInputDialog.getDouble(
                None, "Masa maxima (kDa)", "Masa maxima en kDa (banda mas grande):",
                self.mass_max_kda, mn, 10000000.0, 2,
            )
            if not ok2:
                return
            self.mass_min_kda = mn
            self.mass_max_kda = mx
        else:
            mn, ok1 = QInputDialog.getDouble(
                None, "log(Min) (kDa)", "log10(masa minima en kDa):",
                math.log10(self.mass_min_kda) if self.mass_min_kda > 0 else 0.0,
                -2.0, 10.0, 2,
            )
            if not ok1:
                return
            mx, ok2 = QInputDialog.getDouble(
                None, "log(Max) (kDa)", "log10(masa maxima en kDa):",
                math.log10(self.mass_max_kda) if self.mass_max_kda > 0 else 0.0,
                mn, 10.0, 2,
            )
            if not ok2:
                return
            self.mass_min_kda = 10.0 ** mn
            self.mass_max_kda = 10.0 ** mx
        self._build_scale()
        for lane in self.lane_items:
            lane.update_labels()

    def _toggle_labels(self):
        self.show_labels = not self.show_labels
        for lane in self.lane_items:
            lane.update_labels()

    def _toggle_scale(self):
        self.show_scale = not self.show_scale
        self._build_scale()

    def to_dict(self) -> dict:
        lanes_data = []
        for lane in self.lane_items:
            lanes_data.append(lane.to_dict())
        return {
            "type": "GelElectrophoresisItem",
            "pos": (self.pos().x(), self.pos().y()),
            "lanes": self.num_lanes,
            "width": self.gel_width,
            "height": self.gel_height,
            "show_labels": self.show_labels,
            "show_scale": self.show_scale,
            "scale_unit": self.scale_unit,
            "scale_min": self.scale_min,
            "scale_max": self.scale_max,
            "mass_min_kda": self.mass_min_kda,
            "mass_max_kda": self.mass_max_kda,
            "lanes_data": lanes_data,
        }

    def load_dict(self, data: dict, scene=None):
        self.setPos(data.get("pos", (0, 0))[0], data.get("pos", (0, 0))[1])
        self.show_labels = data.get("show_labels", True)
        self.show_scale = data.get("show_scale", True)
        self.scale_unit = data.get("scale_unit", "Distance")
        self.scale_min = data.get("scale_min", 0.0)
        self.scale_max = data.get("scale_max", 1.0)
        self.mass_min_kda = data.get("mass_min_kda", 10.0)
        self.mass_max_kda = data.get("mass_max_kda", 250.0)
        self.num_lanes = data.get("lanes", 5)
        self.gel_width = data.get("width", 280.0)
        self.gel_height = data.get("height", 320.0)

        for lane in self.lane_items:
            for band, lbl in lane.bands:
                if self.scene():
                    self.scene().removeItem(band)
                    self.scene().removeItem(lbl)
            if self.scene():
                self.scene().removeItem(lane)
        self.lane_items.clear()

        for line in self._scale_lines:
            if line.scene():
                line.scene().removeItem(line)
        for label in self._scale_labels:
            if label.scene():
                label.scene().removeItem(label)
        self._scale_lines.clear()
        self._scale_labels.clear()

        lane_w = self.gel_width / self.num_lanes
        well_height = 14.0
        well_width = lane_w * 0.65

        path = QPainterPath()
        path.moveTo(0, 0)
        for i in range(self.num_lanes):
            cx = lane_w * i + lane_w / 2
            path.lineTo(cx - well_width / 2, 0)
            path.lineTo(cx - well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, well_height)
            path.lineTo(cx + well_width / 2, 0)
        path.lineTo(self.gel_width, 0)
        path.lineTo(self.gel_width, self.gel_height)
        path.lineTo(0, self.gel_height)
        path.closeSubpath()

        self.outline.setRect(0, 0, self.gel_width, self.gel_height)
        self.wells_path.setPath(path)

        lanes_data = data.get("lanes_data", [])
        for i in range(self.num_lanes):
            lane = GelLaneItem(self, index=i, width=lane_w, height=self.gel_height)
            lane.setPos(lane_w * i + lane_w / 2, 0)
            if i < len(lanes_data):
                lane.load_dict(lanes_data[i], scene=scene)
            else:
                lane.add_band(dist=0.0, scene=scene)
            self.lane_items.append(lane)
        self._build_scale()

    def to_json(self) -> dict:
        return self.to_dict()

    @classmethod
    def from_json(cls, data: dict) -> "GelElectrophoresisItem":
        item = cls(
            lanes=data.get("lanes", 5),
            width=data.get("width", 280.0),
            height=data.get("height", 320.0),
        )
        item.load_dict(data)
        return item


# =============================================================================
# Factory
# =============================================================================

class PlateItem:
    """Factory para deserializar placas TLC y geles."""

    @staticmethod
    def from_json(data: dict):
        t = data.get("type")
        if t == "TLCPlateItem":
            return TLCPlateItem.from_json(data)
        elif t == "GelElectrophoresisItem":
            return GelElectrophoresisItem.from_json(data)
        return None
