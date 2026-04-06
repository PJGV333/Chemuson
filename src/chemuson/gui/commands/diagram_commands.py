from __future__ import annotations

from PyQt6.QtCore import QRectF
from PyQt6.QtGui import QUndoCommand

from chemuson.gui.diagram_models import SemanticDiagram

class AddEnergyDiagramItemCommand(QUndoCommand):
    """Comando para añadir un diagrama de energia persistente."""

    def __init__(self, view, item) -> None:
        super().__init__("Add energy diagram")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_energy_diagram_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_energy_diagram_item(self._item)

class TransformEnergyDiagramItemsCommand(QUndoCommand):
    """Comando para mover/escalar diagramas de energia."""

    def __init__(self, view, before: dict, after: dict, text: str = "Transform energy diagrams") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        pos, width, height, rotation = snapshot
        item.set_display_rect(QRectF(float(pos.x()), float(pos.y()), float(width), float(height)))
        item.setRotation(float(rotation))

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class ConfigureEnergyDiagramItemsCommand(QUndoCommand):
    """Comando para cambiar configuracion visual/logica de diagramas de energia."""

    def __init__(self, view, before: dict, after: dict, text: str = "Configure energy diagrams") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        item.apply_config_payload(snapshot)

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class AddCompositeDiagramItemCommand(QUndoCommand):
    """Comando para añadir un diagrama semántico compuesto."""

    def __init__(self, view, item) -> None:
        super().__init__("Add semantic diagram")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_semantic_diagram_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_semantic_diagram_item(self._item)

class EditSemanticDiagramCommand(QUndoCommand):
    """Comando para reconstruir un diagrama semántico existente."""

    def __init__(
        self,
        view,
        item,
        before_payload: dict,
        after_payload: dict,
        text: str = "Edit semantic diagram",
    ) -> None:
        super().__init__(text)
        self._view = view
        self._item = item
        self._before_payload = dict(before_payload or {})
        self._after_payload = dict(after_payload or {})

    def _apply_payload(self, payload: dict) -> None:
        if self._item is None or self._item.scene() is not self._view.scene:
            return
        diagram_payload = dict(payload.get("semantic_diagram", {}) or {})
        self._item.rebuild_from_semantic_diagram(
            SemanticDiagram.from_json_dict(diagram_payload)
        )
        if {"x", "y", "width", "height"} <= set(payload):
            self._item.set_display_rect(
                QRectF(
                    float(payload.get("x", 0.0)),
                    float(payload.get("y", 0.0)),
                    float(payload.get("width", 1.0)),
                    float(payload.get("height", 1.0)),
                )
            )
        if "rotation" in payload:
            self._item.setRotation(float(payload.get("rotation", 0.0)))
        if "z" in payload:
            self._item.setZValue(float(payload.get("z", self._item.zValue())))
        if "opacity" in payload:
            setter = getattr(self._view, "set_graphics_item_opacity", None)
            if callable(setter):
                setter(self._item, payload.get("opacity"))
            else:
                raw = payload.get("opacity")
                self._item.setOpacity(1.0 if raw is None else float(raw))
        self._item.setSelected(True)
        sync_selection = getattr(self._view, "_sync_selection_from_scene", None)
        if callable(sync_selection):
            sync_selection()
        update_overlay = getattr(self._view, "_update_selection_overlay", None)
        if callable(update_overlay):
            update_overlay()

    def redo(self) -> None:
        self._apply_payload(self._after_payload)

    def undo(self) -> None:
        self._apply_payload(self._before_payload)

class AddOrbitalItemCommand(QUndoCommand):
    """Comando para añadir un orbital vectorial persistente."""

    def __init__(self, view, item) -> None:
        super().__init__("Add orbital")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_orbital_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_orbital_item(self._item)

class TransformOrbitalItemsCommand(QUndoCommand):
    """Comando para mover/escalar/rotar orbitales por anchors."""

    def __init__(self, view, before: dict, after: dict, text: str = "Transform orbitals") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        anchor0, anchor1 = snapshot
        item.set_anchors(anchor0, anchor1)

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class StyleOrbitalItemsCommand(QUndoCommand):
    """Comando para aplicar overrides visuales por lóbulo en orbitales."""

    def __init__(self, view, before: dict, after: dict, text: str = "Style orbitals") -> None:
        super().__init__(text)
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        item.set_part_styles(snapshot)

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()
