from __future__ import annotations

from PyQt6.QtGui import QFont, QUndoCommand

class MoveTextItemsCommand(QUndoCommand):
    """Comando para mover elementos de texto en el lienzo."""

    def __init__(self, view, before: dict, after: dict):
        """Inicializa el comando de movimiento de textos."""
        super().__init__("Move text items")
        self._view = view
        self._before = before # {item: (pos, rotation)}
        self._after = after   # {item: (pos, rotation)}

    def redo(self):
        """Aplica nuevas posiciones/rotaciones de texto."""
        for item, (pos, rot) in self._after.items():
            item.setPos(pos)
            item.setRotation(rot)
        self._view._update_selection_overlay()

    def undo(self):
        """Revierte a las posiciones/rotaciones anteriores."""
        for item, (pos, rot) in self._before.items():
            item.setPos(pos)
            item.setRotation(rot)
        self._view._update_selection_overlay()

class ScaleTextItemsCommand(QUndoCommand):
    """Comando para escalar textos libres conservando formato."""

    def __init__(self, view, before: dict, after: dict) -> None:
        super().__init__("Scale text items")
        self._view = view
        self._before = before
        self._after = after

    @staticmethod
    def _apply_snapshot(item, snapshot) -> None:
        pos, rot, font_str, text_width = snapshot
        font = QFont()
        if font_str:
            font.fromString(font_str)
            item.setFont(font)
        item.setPos(pos)
        item.setRotation(rot)
        item.setTextWidth(float(text_width))

    def redo(self) -> None:
        for item, snapshot in self._after.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

    def undo(self) -> None:
        for item, snapshot in self._before.items():
            self._apply_snapshot(item, snapshot)
        self._view._update_selection_overlay()

class FormatTextItemsCommand(QUndoCommand):
    """Comando para aplicar formato tipográfico con undo/redo."""

    def __init__(
        self,
        view,
        items: list,
        settings_before: dict,
        settings_after: dict,
        property_name: str,
    ) -> None:
        super().__init__("Format text")
        self._view = view
        self._items = [item for item in items if item is not None and item.scene() is view.scene]
        self._settings_before = dict(settings_before)
        self._settings_after = dict(settings_after)
        self._property_name = property_name
        self._before = {
            item: view._text_item_format_snapshot(item)
            for item in self._items
        }
        self._after = None

    def redo(self) -> None:
        self._view._current_text_settings.update(self._settings_after)
        if self._after is None:
            for item in self._items:
                self._view._apply_text_settings(item, self._property_name)
            self._after = {
                item: self._view._text_item_format_snapshot(item)
                for item in self._items
            }
        else:
            for item, snapshot in self._after.items():
                self._view._restore_text_item_format_snapshot(item, snapshot)
        self._view.sync_text_selection_state()

    def undo(self) -> None:
        self._view._current_text_settings.update(self._settings_before)
        for item, snapshot in self._before.items():
            self._view._restore_text_item_format_snapshot(item, snapshot)
        self._view.sync_text_selection_state()

class AddTextItemCommand(QUndoCommand):
    """Comando para añadir un elemento de texto."""

    def __init__(self, view, item) -> None:
        """Inicializa el comando de texto."""
        super().__init__("Add text")
        self._view = view
        self._item = item

    def redo(self) -> None:
        """Reintroduce el texto en el lienzo."""
        if self._item is not None:
            self._view.readd_text_item(self._item)

    def undo(self) -> None:
        """Elimina el texto añadido."""
        if self._item is not None:
            self._view.remove_text_item(self._item)
