from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QUndoCommand

class ChangeArrowStrokeCommand(QUndoCommand):
    """Comando para cambiar el grosor de una flecha seleccionada."""

    def __init__(self, view, item, new_stroke_px: Optional[float]) -> None:
        super().__init__("Change arrow thickness")
        self._view = view
        self._item = item
        self._new_stroke = new_stroke_px
        self._old_stroke = item.stroke_px() if hasattr(item, "stroke_px") else None

    def _apply(self, stroke_px: Optional[float]) -> None:
        if self._item is None or not hasattr(self._item, "set_stroke_px"):
            return
        self._item.set_stroke_px(stroke_px)
        if hasattr(self._view, "_update_selection_overlay"):
            self._view._update_selection_overlay()

    def redo(self) -> None:
        self._apply(self._new_stroke)

    def undo(self) -> None:
        self._apply(self._old_stroke)

class ChangeBracketStrokeCommand(QUndoCommand):
    """Comando para cambiar el grosor de un corchete/llave/paréntesis."""

    def __init__(self, view, item, new_stroke_px: Optional[float]) -> None:
        super().__init__("Change bracket thickness")
        self._view = view
        self._item = item
        self._new_stroke = new_stroke_px
        self._old_stroke = item.stroke_px() if hasattr(item, "stroke_px") else None

    def _apply(self, stroke_px: Optional[float]) -> None:
        if self._item is None or not hasattr(self._item, "set_stroke_px"):
            return
        self._item.set_stroke_px(stroke_px)
        if hasattr(self._view, "_update_selection_overlay"):
            self._view._update_selection_overlay()

    def redo(self) -> None:
        self._apply(self._new_stroke)

    def undo(self) -> None:
        self._apply(self._old_stroke)

class AddArrowCommand(QUndoCommand):
    """Comando para añadir una flecha de anotación."""

    def __init__(
        self,
        view,
        start: QPointF,
        end: QPointF,
        kind: str,
        curve_factor: float | None = None,
        stroke_px: float | None = None,
        opacity: Optional[float] = None,
    ) -> None:
        """Inicializa el comando de flecha."""
        super().__init__("Add arrow")
        self._view = view
        self._start = QPointF(start)
        self._end = QPointF(end)
        self._kind = kind
        self._curve_factor = curve_factor
        self._stroke_px = stroke_px
        self._opacity = opacity
        self._item = None

    def redo(self) -> None:
        """Crea o reintroduce la flecha."""
        if self._item is None:
            self._item = self._view.add_arrow_item(
                self._start,
                self._end,
                self._kind,
                curve_factor=self._curve_factor,
                stroke_px=self._stroke_px,
                opacity=self._opacity,
            )
        else:
            self._view.readd_arrow_item(
                self._item,
                self._start,
                self._end,
                self._kind,
                curve_factor=self._curve_factor,
                opacity=self._opacity,
            )

    def undo(self) -> None:
        """Elimina la flecha añadida."""
        if self._item is not None:
            self._view.remove_arrow_item(self._item)

    @property
    def item(self):
        """Devuelve el item de flecha creado por el comando."""
        return self._item

class AddBracketCommand(QUndoCommand):
    """Comando para añadir corchetes/llaves de anotación."""

    def __init__(
        self,
        view,
        rect: QRectF,
        kind: str,
        padding: float | None = None,
        stroke_px: float | None = None,
        opacity: Optional[float] = None,
    ) -> None:
        """Inicializa el comando de corchetes."""
        super().__init__("Add brackets")
        self._view = view
        self._rect = QRectF(rect)
        self._kind = kind
        self._padding = 8.0 if padding is None else float(padding)
        self._stroke_px = stroke_px
        self._opacity = opacity
        self._item = None

    def redo(self) -> None:
        """Crea o reintroduce los corchetes."""
        if self._item is None:
            self._item = self._view.add_bracket_item(
                self._rect,
                self._kind,
                padding=self._padding,
                stroke_px=self._stroke_px,
                opacity=self._opacity,
            )
        else:
            self._view.readd_bracket_item(
                self._item,
                self._rect,
                self._kind,
                padding=self._padding,
                stroke_px=self._stroke_px,
                opacity=self._opacity,
            )

    def undo(self) -> None:
        """Elimina los corchetes añadidos."""
        if self._item is not None:
            self._view.remove_bracket_item(self._item)

    @property
    def item(self):
        """Devuelve el item de corchete creado por el comando."""
        return self._item

class AddImageItemCommand(QUndoCommand):
    """Comando para añadir una imagen anotada."""

    def __init__(self, view, item) -> None:
        super().__init__("Add image")
        self._view = view
        self._item = item

    def redo(self) -> None:
        if self._item is not None:
            self._view.readd_image_item(self._item)

    def undo(self) -> None:
        if self._item is not None:
            self._view.remove_image_item(self._item)

class AddWavyAnchorCommand(QUndoCommand):
    """Comando para añadir un ancla ondulada."""

    def __init__(self, view, item) -> None:
        """Inicializa el comando de ancla ondulada."""
        super().__init__("Add wavy anchor")
        self._view = view
        self._item = item

    def redo(self) -> None:
        """Reintroduce el ancla ondulada."""
        if self._item is not None:
            self._view.readd_wavy_anchor_item(self._item)

    def undo(self) -> None:
        """Elimina el ancla ondulada añadida."""
        if self._item is not None:
            self._view.remove_wavy_anchor_item(self._item)
