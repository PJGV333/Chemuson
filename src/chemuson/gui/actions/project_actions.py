"""Creación de acciones de proyecto/estructura."""

from PyQt6.QtGui import QAction, QActionGroup, QKeySequence

from chemuson.gui.canvas import BRANCH_ROTATION_STEP_DEG, FRAGMENT_ROTATION_STEP_DEG


def create_project_actions(window) -> None:
    """Inicializa acciones de vista/estructura/proyecto en ``window``."""
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

    window.action_numbering_enabled = QAction("Mostrar numeración", window)
    window.action_numbering_enabled.setCheckable(True)
    window.action_numbering_enabled.setChecked(False)
    window.action_numbering_enabled.triggered.connect(window._on_toggle_numbering)

    window.action_numbering_mode_atoms = QAction("Numerar átomos", window)
    window.action_numbering_mode_atoms.setCheckable(True)
    window.action_numbering_mode_atoms.triggered.connect(
        lambda checked=False: window._on_set_numbering_mode("atoms")
    )
    window.action_numbering_mode_structures = QAction("Numerar estructuras", window)
    window.action_numbering_mode_structures.setCheckable(True)
    window.action_numbering_mode_structures.triggered.connect(
        lambda checked=False: window._on_set_numbering_mode("structures")
    )
    window.action_numbering_mode_both = QAction("Numerar ambos", window)
    window.action_numbering_mode_both.setCheckable(True)
    window.action_numbering_mode_both.triggered.connect(
        lambda checked=False: window._on_set_numbering_mode("both")
    )
    window._numbering_mode_group = QActionGroup(window)
    window._numbering_mode_group.setExclusive(True)
    window._numbering_mode_group.addAction(window.action_numbering_mode_atoms)
    window._numbering_mode_group.addAction(window.action_numbering_mode_structures)
    window._numbering_mode_group.addAction(window.action_numbering_mode_both)
    window.action_numbering_mode_atoms.setChecked(True)

    window.action_numbering_recalculate = QAction("Recalcular numeración", window)
    window.action_numbering_recalculate.triggered.connect(window._on_recalculate_numbering)

    window.action_numbering_export = QAction("Incluir numeración en exportación", window)
    window.action_numbering_export.setCheckable(True)
    window.action_numbering_export.setChecked(True)
    window.action_numbering_export.triggered.connect(window._on_toggle_numbering_export)

    window.action_clean_2d = QAction("Limpiar 2D", window)
    window.action_clean_2d.triggered.connect(window._on_clean_2d)
    window.action_clean_2d_full = QAction("Limpiar 2D (1 paso)", window)
    window.action_clean_2d_full.setShortcut(QKeySequence("Ctrl+K"))
    window.action_clean_2d_full.triggered.connect(window._on_clean_2d_full)

    window.action_rotate_left = QAction("Girar 90° a la izquierda", window)
    window.action_rotate_left.triggered.connect(lambda: window._on_rotate_selection(-90.0))

    window.action_rotate_right = QAction("Girar 90° a la derecha", window)
    window.action_rotate_right.triggered.connect(lambda: window._on_rotate_selection(90.0))

    window.action_flip_horizontal = QAction("Giro 180° horizontal", window)
    window.action_flip_horizontal.triggered.connect(window._on_flip_horizontal)

    window.action_flip_vertical = QAction("Giro 180° vertical", window)
    window.action_flip_vertical.triggered.connect(window._on_flip_vertical)

    window.action_branch_rotate_minus = QAction(f"Girar rama -{int(BRANCH_ROTATION_STEP_DEG)}°", window)
    window.action_branch_rotate_minus.setShortcut(QKeySequence("Ctrl+Alt+Left"))
    window.action_branch_rotate_minus.triggered.connect(
        lambda: window._on_rotate_branch(-BRANCH_ROTATION_STEP_DEG)
    )

    window.action_branch_rotate_plus = QAction(f"Girar rama +{int(BRANCH_ROTATION_STEP_DEG)}°", window)
    window.action_branch_rotate_plus.setShortcut(QKeySequence("Ctrl+Alt+Right"))
    window.action_branch_rotate_plus.triggered.connect(
        lambda: window._on_rotate_branch(BRANCH_ROTATION_STEP_DEG)
    )

    window.action_branch_invert = QAction("Invertir rama (180°)", window)
    window.action_branch_invert.setShortcut(QKeySequence("Ctrl+Alt+I"))
    window.action_branch_invert.triggered.connect(window._on_invert_branch)

    window.action_branch_auto_arrange = QAction("Autoacomodar rama", window)
    window.action_branch_auto_arrange.setShortcut(QKeySequence("Ctrl+Alt+A"))
    window.action_branch_auto_arrange.triggered.connect(window._on_auto_arrange_branch)

    window.action_fragment_pivot_set = QAction("Definir átomo pivote desde selección", window)
    window.action_fragment_pivot_set.triggered.connect(window._on_set_fragment_pivot)

    window.action_fragment_pivot_clear = QAction("Limpiar átomo pivote", window)
    window.action_fragment_pivot_clear.triggered.connect(window._on_clear_fragment_pivot)

    window.action_fragment_rotate_minus = QAction(f"Girar fragmento -{int(FRAGMENT_ROTATION_STEP_DEG)}°", window)
    window.action_fragment_rotate_minus.triggered.connect(
        lambda: window._on_rotate_fragment(-FRAGMENT_ROTATION_STEP_DEG)
    )

    window.action_fragment_rotate_plus = QAction(f"Girar fragmento +{int(FRAGMENT_ROTATION_STEP_DEG)}°", window)
    window.action_fragment_rotate_plus.triggered.connect(
        lambda: window._on_rotate_fragment(FRAGMENT_ROTATION_STEP_DEG)
    )

    window.action_fragment_invert = QAction("Invertir fragmento (180°)", window)
    window.action_fragment_invert.triggered.connect(window._on_invert_fragment)
