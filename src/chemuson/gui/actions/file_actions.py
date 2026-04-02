"""Creación de acciones relacionadas con gestión de archivos."""

from PyQt6.QtGui import QAction, QKeySequence


def create_file_actions(window) -> None:
    """Inicializa acciones de archivo en ``window``."""
    window.action_new = QAction("Nuevo", window)
    window.action_new.setShortcut(QKeySequence.StandardKey.New)
    window.action_new.triggered.connect(window._on_file_new)

    window.action_open = QAction("Abrir...", window)
    window.action_open.setShortcut(QKeySequence.StandardKey.Open)
    window.action_open.triggered.connect(window._on_file_open)

    window.action_save = QAction("Guardar", window)
    window.action_save.setShortcuts([
        QKeySequence.StandardKey.Save,
        QKeySequence("Ctrl+G"),
    ])
    window.action_save.triggered.connect(window._on_file_save)

    window.action_recovery_center = QAction("Centro de recuperación...", window)
    window.action_recovery_center.triggered.connect(window._on_open_recovery_center)

    window.action_quit = QAction("Salir", window)
    window.action_quit.setShortcut(QKeySequence.StandardKey.Quit)
    window.action_quit.triggered.connect(window.close)
