"""Acciones de vista/proyección del documento."""

from PyQt6.QtGui import QAction, QKeySequence


def create_view_actions(window) -> None:
    window.action_show_carbons = QAction("Mostrar carbonos", window)
    window.action_show_carbons.setCheckable(True)
    window.action_show_carbons.setChecked(False)
    window.action_show_carbons.triggered.connect(window._on_toggle_carbons)

    window.action_show_hydrogens = QAction("Mostrar hidrógenos", window)
    window.action_show_hydrogens.setCheckable(True)
    window.action_show_hydrogens.setChecked(False)
    window.action_show_hydrogens.triggered.connect(window._on_toggle_hydrogens)

    window.action_aromatic_circles = QAction("Aromáticos como círculos", window)
    window.action_aromatic_circles.setCheckable(True)
    window.action_aromatic_circles.setChecked(False)
    window.action_aromatic_circles.triggered.connect(window._on_toggle_aromatic_circles)

    window.action_zoom_in = QAction("Zoom +", window)
    window.action_zoom_in.setShortcut(QKeySequence.StandardKey.ZoomIn)
    window.action_zoom_in.triggered.connect(window._on_zoom_in)

    window.action_zoom_out = QAction("Zoom -", window)
    window.action_zoom_out.setShortcut(QKeySequence.StandardKey.ZoomOut)
    window.action_zoom_out.triggered.connect(window._on_zoom_out)

    window.action_zoom_reset = QAction("Zoom 100%", window)
    window.action_zoom_reset.setShortcut("Ctrl+0")
    window.action_zoom_reset.triggered.connect(window._on_zoom_reset)

    window.action_show_main_toolbar_aux = QAction(
        "Mostrar copiar/pegar/zoom en barra superior",
        window,
    )
    window.action_show_main_toolbar_aux.setCheckable(True)
    window.action_show_main_toolbar_aux.setChecked(False)
    window.action_show_main_toolbar_aux.toggled.connect(window._on_toggle_main_toolbar_aux)

    window.action_rules = QAction("Reglas", window)
    window.action_rules.setCheckable(True)
    window.action_rules.setChecked(False)
    window.action_rules.triggered.connect(window._on_toggle_rules)

    window.action_crosshair = QAction("Cuadrícula", window)
    window.action_crosshair.setCheckable(True)
    window.action_crosshair.setChecked(False)
    window.action_crosshair.triggered.connect(window._on_toggle_crosshair)
