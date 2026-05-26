"""Acciones de estructura/transformación molecular."""

from PyQt6.QtGui import QAction, QKeySequence

from chemuson.gui.canvas import BRANCH_ROTATION_STEP_DEG, FRAGMENT_ROTATION_STEP_DEG


def create_structure_actions(window) -> None:
    window.action_clean_2d = QAction("Limpiar 2D", window)
    window.action_clean_2d.triggered.connect(window._on_clean_2d)
    window.action_clean_2d_full = QAction("Limpiar 2D (1 paso)", window)
    window.action_clean_2d_full.setShortcut(QKeySequence("Ctrl+K"))
    window.action_clean_2d_full.triggered.connect(window._on_clean_2d_full)
    window.action_clean_2d_publication = QAction("Limpiar 2D para publicación", window)
    window.action_clean_2d_publication.setShortcut(QKeySequence("Ctrl+Shift+K"))
    window.action_clean_2d_publication.triggered.connect(window._on_clean_2d_publication)

    window.action_validate_structure = QAction("Validar valencias", window)
    window.action_validate_structure.setShortcut(QKeySequence("Ctrl+Shift+V"))
    window.action_validate_structure.triggered.connect(window._on_validate_structure)

    window.action_validation_next = QAction("Siguiente error de valencia", window)
    window.action_validation_next.setShortcut(QKeySequence("F8"))
    window.action_validation_next.triggered.connect(lambda: window._on_navigate_validation_issue(1))

    window.action_validation_previous = QAction("Error de valencia anterior", window)
    window.action_validation_previous.setShortcut(QKeySequence("Shift+F8"))
    window.action_validation_previous.triggered.connect(lambda: window._on_navigate_validation_issue(-1))

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
