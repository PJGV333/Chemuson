"""Creación de acciones relacionadas con actualizaciones."""

from PyQt6.QtGui import QAction


def create_update_actions(window) -> None:
    """Inicializa acciones de actualización en ``window``."""
    window.action_check_updates_now = QAction("Buscar actualizaciones...", window)
    window.action_check_updates_now.triggered.connect(window._on_check_updates_now)
