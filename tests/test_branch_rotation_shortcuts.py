"""Shortcuts de ventana para rotación de rama: funcionan sin abrir el menú.

Los QActions históricos "Girar rama ±60° / Invertir rama / Autoacomodar rama"
(``Editar -> Rotar``) vivían únicamente asociados al QMenu de la barra de
menús (oculta en la UI moderna, ``menuBar().setVisible(False)``). Con la
barra oculta, el shortcut map de esos QMenus NO está activo cuando el foco
está en el lienzo (ni en ningún otro hijo de la ventana), por lo que
``Ctrl+Alt+Left`` / ``Ctrl+Alt+Right`` / ``Ctrl+Alt+I`` / ``Ctrl+Alt+A``
no se activaban con el menú cerrado; el evento caía en el nudge de flechas
del canvas (traslación de 1 px, imperceptible) o en el vacío.

El arreglo registra la misma QAction en la ventana
(``Qt.ShortcutContext.WindowShortcut`` + ``window.addAction``), sin duplicar
acciones ni conexiones ``triggered`` (mismo patrón de atajo de ventana que las
acciones de la barra superior y la paleta de comandos).
"""

from __future__ import annotations

import contextlib
import math

import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QAction
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QLineEdit

from chemuson.gui.geom import angle_deg, endpoint_from_angle_len

_CTRL_ALT = Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.AltModifier


@contextlib.contextmanager
def _window():
    """Ventana real con foco en el lienzo (estado natural de uso); garantiza
    limpieza del undo stack antes de cerrar (sin diálogo de descarte)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    win.canvas.setFocus()
    QApplication.processEvents()
    try:
        yield win
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def _add_chain_with_pivot_bond(win):
    """Molécula acíclica A-B-C-(CD)-D-E con exactamente el enlace CD
    acíclico seleccionado (misma geometría de ``test_branch_reorientation``)."""
    canvas = win.canvas
    atoms = [
        canvas.model.add_atom("C", x, y)
        for x, y in [(80.0, 100.0), (120.0, 100.0), (160.0, 100.0), (200.0, 60.0), (240.0, 60.0)]
    ]
    ids = [atom.id for atom in atoms]
    canvas.model.add_bond(ids[0], ids[1], order=1)
    canvas.model.add_bond(ids[1], ids[2], order=1)
    bond_cd = canvas.model.add_bond(ids[2], ids[3], order=1)
    canvas.model.add_bond(ids[3], ids[4], order=1)
    canvas._rebuild_items_from_model()
    canvas.bond_items[bond_cd.id].setSelected(True)
    canvas._sync_selection_from_scene()
    return ids


def _add_obstacle_branch(win):
    """Geometría de autoarrange con obstáculo: el lado en movimiento
    (N + tail) debe ser la rama MENOR del enlace seleccionado (el contexto
    por defecto autoacomoda la rama menor), y el obstáculo debe hacer que
    la orientación actual (120°) sea peor que la alternativa (240°).
    El enlace pivote center-N queda seleccionado."""
    canvas = win.canvas
    center = canvas.model.add_atom("C", 200.0, 200.0)
    carbonyl_o = canvas.model.add_atom("O", 240.0, 200.0)
    extra_sub = canvas.model.add_atom("Cl", 160.0, 200.0)
    moving_pos = endpoint_from_angle_len(QPointF(200.0, 200.0), 120.0, 40.0)
    moving = canvas.model.add_atom("N", moving_pos.x(), moving_pos.y())
    tail_pos = endpoint_from_angle_len(moving_pos, 120.0, 40.0)
    tail = canvas.model.add_atom("C", tail_pos.x(), tail_pos.y())
    obstacle = canvas.model.add_atom("Cl", tail_pos.x() + 2.0, tail_pos.y() + 1.0)

    canvas.model.add_bond(center.id, carbonyl_o.id, order=2)
    canvas.model.add_bond(center.id, extra_sub.id, order=1)
    pivot = canvas.model.add_bond(center.id, moving.id, order=1)
    canvas.model.add_bond(moving.id, tail.id, order=1)
    canvas._rebuild_items_from_model()
    canvas.bond_items[pivot.id].setSelected(True)
    canvas._sync_selection_from_scene()
    return center.id, moving.id, tail.id, obstacle.id


def _coords(canvas, atom_ids):
    return {
        atom_id: (canvas.model.get_atom(atom_id).x, canvas.model.get_atom(atom_id).y)
        for atom_id in atom_ids
    }


def _max_displacement(before, after) -> float:
    return max(
        math.hypot(after[i][0] - before[i][0], after[i][1] - before[i][1])
        for i in before
    )


def _send_key(win, widget, key: Qt.Key) -> None:
    """Envía Ctrl+Alt+<key> al widget (con foco), sin abrir ningún menú."""
    QTest.keyClick(widget, key, _CTRL_ALT)
    QApplication.processEvents()


def test_branch_rotation_shortcuts_registered_on_window():
    """Wiring: las cuatro acciones de rama registran el shortcut en la
    ventana (mismo patrón que Clean2D) y el mismo QAction histórico sigue
    viviendo en el menú Rotar, sin QActions duplicados."""
    with _window() as win:
        actions = (
            win.action_branch_rotate_minus,
            win.action_branch_rotate_plus,
            win.action_branch_invert,
            win.action_branch_auto_arrange,
        )
        for action in actions:
            assert action.shortcutContext() == Qt.ShortcutContext.WindowShortcut
            # ``window.addAction`` deja la acción en la lista de la ventana
            # (el discriminador real: solo asociada al menú no aparece aquí).
            assert action in win.actions()
        # El mismo QAction histórico sigue viviendo en el menú Rotar.
        rotate_menu = None
        for action in win.menuBar().actions():
            submenu = action.menu()
            if submenu is not None and submenu.title() == "Editar":
                for subaction in submenu.actions():
                    if subaction.menu() is not None and subaction.menu().title() == "Rotar":
                        rotate_menu = subaction.menu()
        assert rotate_menu is not None
        for action in actions:
            assert action in rotate_menu.actions()
        # Sin duplicados: los 4 objetos son los mismos QActions históricos
        # (mismas identidades) y cada QAction de la ventana es único.
        assert len({id(a) for a in actions}) == 4
        all_actions = win.findChildren(QAction)
        assert len(all_actions) == len({id(a) for a in all_actions})


def test_ctrl_alt_right_rotates_branch_with_menu_closed():
    """E2E: Ctrl+Alt+Right con el menú cerrado rota la rama ±60° (no el
    nudge de 1 px del canvas) y soporta undo/redo: una sola
    transformación por pulsación (sin doble ejecución)."""
    with _window() as win:
        atom_ids = _add_chain_with_pivot_bond(win)
        canvas = win.canvas
        before = _coords(canvas, atom_ids)
        _send_key(win, canvas, Qt.Key.Key_Right)
        after = _coords(canvas, atom_ids)
        assert after != before
        # Una rotación de 60° desplaza átomos muchos px; el nudge del canvas
        # solo movería 1 px.
        assert _max_displacement(before, after) > 5.0
        # Una sola ejecución (sin doble transformación del shortcut).
        assert canvas.undo_stack.count() == 1
        canvas.undo_stack.undo()
        restored = _coords(canvas, atom_ids)
        for atom_id in atom_ids:
            assert restored[atom_id] == pytest.approx(before[atom_id])
        canvas.undo_stack.redo()
        redone = _coords(canvas, atom_ids)
        for atom_id in atom_ids:
            assert redone[atom_id] == pytest.approx(after[atom_id])


def test_ctrl_alt_left_rotates_branch_with_menu_closed():
    """E2E: Ctrl+Alt+Left con el menú cerrado también rota la rama."""
    with _window() as win:
        atom_ids = _add_chain_with_pivot_bond(win)
        canvas = win.canvas
        before = _coords(canvas, atom_ids)
        _send_key(win, canvas, Qt.Key.Key_Left)
        after = _coords(canvas, atom_ids)
        assert after != before
        assert _max_displacement(before, after) > 5.0
        assert canvas.undo_stack.count() == 1


def test_ctrl_alt_i_inverts_branch_with_menu_closed():
    """E2E: Ctrl+Alt+I invierte la rama 180° (geometría esperada conocida
    de ``test_invert_selected_branch_uses_smaller_side_and_preserves_lengths``)."""
    with _window() as win:
        atom_ids = _add_chain_with_pivot_bond(win)
        canvas = win.canvas
        _send_key(win, canvas, Qt.Key.Key_I)
        d = canvas.model.get_atom(atom_ids[3])
        e = canvas.model.get_atom(atom_ids[4])
        assert (d.x, d.y) == pytest.approx((120.0, 140.0))
        assert (e.x, e.y) == pytest.approx((80.0, 140.0))
        assert canvas.undo_stack.count() == 1


def test_ctrl_alt_a_auto_arranges_branch_with_menu_closed():
    """E2E: Ctrl+Alt+A autoacomoda la rama hacia el lado menos congestionado
    (el ángulo del átomo en movimiento pasa de ~120° a ~240°)."""
    with _window() as win:
        center_id, moving_id, _tail_id, _obstacle_id = _add_obstacle_branch(win)
        canvas = win.canvas

        def _angle() -> float:
            c = canvas.model.get_atom(center_id)
            m = canvas.model.get_atom(moving_id)
            return angle_deg(QPointF(c.x, c.y), QPointF(m.x, m.y))

        assert _angle() == pytest.approx(120.0, abs=0.2)
        _send_key(win, canvas, Qt.Key.Key_A)
        assert _angle() == pytest.approx(240.0, abs=0.4)
        assert canvas.undo_stack.count() == 1


def test_menu_path_still_works_after_fix():
    """La ruta de menú sigue funcionando: el mismo QAction histórico (el que
    vive en ``Editar -> Rotar``) al activarse rota la rama."""
    with _window() as win:
        atom_ids = _add_chain_with_pivot_bond(win)
        canvas = win.canvas
        before = _coords(canvas, atom_ids)
        win.action_branch_rotate_plus.trigger()
        QApplication.processEvents()
        after = _coords(canvas, atom_ids)
        assert after != before
        assert _max_displacement(before, after) > 5.0
        assert canvas.undo_stack.count() == 1


def test_shortcut_with_non_modal_text_editor_focus_keeps_window_policy():
    """Editor de texto NO modal (hijo de la ventana): la política existente
    de supresión (``_tool_shortcuts_suppressed``) cubre las letras de
    herramienta (que además escribirían en el editor); las combinaciones
    Ctrl+Alt+* no escriben nada y siguen activas por semántica estándar de
    ``WindowShortcut`` (un hijo con foco cuenta como foco de la ventana —
    mismo comportamiento que el ``Ctrl+P`` de la paleta de comandos). El shortcut
    debe activarse aunque el foco esté en el editor, enviando la tecla al
    widget con foco (como haría un pulsado real)."""
    with _window() as win:
        atom_ids = _add_chain_with_pivot_bond(win)
        canvas = win.canvas
        before = _coords(canvas, atom_ids)
        line = QLineEdit(win)
        line.show()
        line.setFocus()
        QApplication.processEvents()
        _send_key(win, line, Qt.Key.Key_Right)
        after = _coords(canvas, atom_ids)
        assert after != before
        line.deleteLater()
