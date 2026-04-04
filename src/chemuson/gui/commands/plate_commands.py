from __future__ import annotations

from ._shared import *

class AddPlateItemCommand(QUndoCommand):
    """Comando para añadir placas de cromatografía y electroforesis."""

    def __init__(self, view, item) -> None:
        super().__init__("Add analysis plate")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_plate_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_plate_item(self._item)

class MovePlateItemsCommand(QUndoCommand):
    """Comando para mover placas."""

    def __init__(self, view, before: dict, after: dict) -> None:
        super().__init__("Move plates")
        self._view = view
        self._before = before
        self._after = after

    def redo(self) -> None:
        for item, (pos, rot) in self._after.items():
            if item.scene() is self._view.scene:
                item.setPos(pos)
                item.setRotation(rot)

    def undo(self) -> None:
        for item, (pos, rot) in self._before.items():
            if item.scene() is self._view.scene:
                item.setPos(pos)
                item.setRotation(rot)

class TransformPlateItemsCommand(QUndoCommand):
    """Comando para escalar y rotar placas."""

    def __init__(self, view, before: dict, after: dict, text="Transform plates") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            if item.scene() is self._view.scene:
                item.set_display_rect(snapshot[0])
                item.setRotation(snapshot[1])

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            if item.scene() is self._view.scene:
                item.set_display_rect(snapshot[0])
                item.setRotation(snapshot[1])

class AddSpotBandCommand(QUndoCommand):
    """Comando para añadir una mancha TLC o banda de gel."""

    def __init__(self, parent_item, spot_or_band, lane_item) -> None:
        super().__init__("Add spot/band")
        self._parent = parent_item
        self._item = spot_or_band
        self._lane = lane_item

    def redo(self) -> None:
        if self._item is not None:
            self._item.setVisible(True)
            if hasattr(self._item, "parentItem") and self._item.parentItem() is None:
                self._item.setParentItem(self._lane)
            if hasattr(self._lane, "rf_labels"):
                if (self._item, None) not in self._lane.rf_labels:
                    from PyQt6.QtGui import QFont, QColor
                    from PyQt6.QtWidgets import QGraphicsTextItem
                    label = QGraphicsTextItem("", self._lane)
                    label.setFont(QFont("Arial", 8))
                    label.setDefaultTextColor(QColor(30, 41, 59))
                    self._lane.rf_labels.append((self._item, label))
                    self._lane.update_rf_labels()
            elif hasattr(self._lane, "bands"):
                if (self._item, None) not in self._lane.bands:
                    from PyQt6.QtGui import QFont, QColor
                    from PyQt6.QtWidgets import QGraphicsTextItem
                    label = QGraphicsTextItem("", self._lane)
                    label.setFont(QFont("Arial", 8))
                    label.setDefaultTextColor(QColor(30, 41, 59))
                    self._lane.bands.append((self._item, label))
                    self._lane.update_labels()

    def undo(self) -> None:
        if self._item is not None:
            scene = self._item.scene()
            if scene:
                scene.removeItem(self._item)
            if hasattr(self._lane, "rf_labels"):
                self._lane.rf_labels = [(s, l) for s, l in self._lane.rf_labels if s is not self._item]
                self._lane.update_rf_labels()
            elif hasattr(self._lane, "bands"):
                self._lane.bands = [(b, l) for b, l in self._lane.bands if b is not self._item]
                self._lane.update_labels()

class RemoveSpotBandCommand(QUndoCommand):
    """Comando para eliminar una mancha TLC o banda de gel."""

    def __init__(self, lane_item, spot_or_band, label_item=None) -> None:
        super().__init__("Remove spot/band")
        self._lane = lane_item
        self._item = spot_or_band
        self._label = label_item
        self._saved_data = None

    def redo(self) -> None:
        if self._item is not None:
            self._saved_data = self._item.to_dict()
            scene = self._item.scene()
            if scene:
                scene.removeItem(self._item)
            if self._label and self._label.scene():
                self._label.scene().removeItem(self._label)
            if hasattr(self._lane, "rf_labels"):
                self._lane.rf_labels = [(s, l) for s, l in self._lane.rf_labels if s is not self._item]
                self._lane.update_rf_labels()
            elif hasattr(self._lane, "bands"):
                self._lane.bands = [(b, l) for b, l in self._lane.bands if b is not self._item]
                self._lane.update_labels()

    def undo(self) -> None:
        if self._saved_data is not None:
            from chemuson.gui.plate_items import TLCSpotItem, GelBandItem
            if isinstance(self._item, TLCSpotItem):
                spot = TLCSpotItem(self._lane)
                spot.load_dict(self._saved_data)
                label = QGraphicsTextItem("", self._lane)
                label.setFont(QFont("Arial", 8))
                label.setDefaultTextColor(QColor(30, 41, 59))
                self._lane.rf_labels.append((spot, label))
                self._lane.update_rf_labels()
            elif isinstance(self._item, GelBandItem):
                band = GelBandItem(self._lane, width=self._saved_data.get("width", 32))
                band.load_dict(self._saved_data)
                label = QGraphicsTextItem("", self._lane)
                label.setFont(QFont("Arial", 8))
                label.setDefaultTextColor(QColor(30, 41, 59))
                self._lane.bands.append((band, label))
                self._lane.update_labels()

class MoveSpotBandCommand(QUndoCommand):
    """Comando para mover una mancha o banda dentro de su carril."""

    def __init__(self, item, before_y: float, after_y: float) -> None:
        super().__init__("Move spot/band")
        self._item = item
        self._before_y = before_y
        self._after_y = after_y

    def redo(self) -> None:
        if self._item is not None:
            self._item.setY(self._after_y)

    def undo(self) -> None:
        if self._item is not None:
            self._item.setY(self._before_y)

class ChangeSpotBandPropertyCommand(QUndoCommand):
    """Comando para cambiar propiedades de mancha/banda (color, opacidad, tamaño)."""

    def __init__(self, item, property_name: str, before_value, after_value) -> None:
        super().__init__(f"Change {property_name}")
        self._item = item
        self._property = property_name
        self._before = before_value
        self._after = after_value

    def redo(self) -> None:
        if self._item is not None:
            self._apply(self._after)

    def undo(self) -> None:
        if self._item is not None:
            self._apply(self._before)

    def _apply(self, value):
        if self._property == "color":
            self._item.set_color(value)
        elif self._property == "opacity":
            self._item.set_opacity(value)
        elif self._property == "size":
            self._item.set_size(*value)

class ChangePlateLabelsCommand(QUndoCommand):
    """Comando para cambiar visibilidad de etiquetas en placa/gel."""

    def __init__(self, plate_item, before: bool, after: bool) -> None:
        super().__init__("Toggle plate labels")
        self._plate = plate_item
        self._before = before
        self._after = after

    def redo(self) -> None:
        self._plate.show_labels = self._after
        for lane in self._plate.lane_items:
            if hasattr(lane, "update_rf_labels"):
                lane.update_rf_labels()
            elif hasattr(lane, "update_labels"):
                lane.update_labels()

    def undo(self) -> None:
        self._plate.show_labels = self._before
        for lane in self._plate.lane_items:
            if hasattr(lane, "update_rf_labels"):
                lane.update_rf_labels()
            elif hasattr(lane, "update_labels"):
                lane.update_labels()
