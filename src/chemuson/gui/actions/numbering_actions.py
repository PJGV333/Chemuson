"""Acciones de numeración."""

from PyQt6.QtGui import QAction, QActionGroup


def create_numbering_actions(window) -> None:
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
