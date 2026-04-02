"""Creación de acciones de edición y portapapeles."""

from PyQt6.QtGui import QAction, QKeySequence


def create_edit_actions(window) -> None:
    """Inicializa acciones de edición en ``window``."""
    window.action_undo = QAction("Deshacer", window)
    window.action_undo.setShortcut(QKeySequence.StandardKey.Undo)

    window.action_redo = QAction("Rehacer", window)
    window.action_redo.setShortcut(QKeySequence.StandardKey.Redo)

    window.action_copy = QAction("Copiar", window)
    window.action_copy.setShortcut(QKeySequence.StandardKey.Copy)

    window.action_cut = QAction("Cortar", window)
    window.action_cut.setShortcut(QKeySequence.StandardKey.Cut)

    window.action_paste = QAction("Pegar", window)
    window.action_paste.setShortcut(QKeySequence.StandardKey.Paste)

    window.action_duplicate = QAction("Duplicar", window)
    window.action_duplicate.setShortcut(QKeySequence("Ctrl+D"))

    window.action_delete = QAction("Eliminar", window)
    window.action_delete.setShortcut(QKeySequence.StandardKey.Delete)

    window.action_edit_electronic_diagram = QAction("Edit Electronic Diagram...", window)
