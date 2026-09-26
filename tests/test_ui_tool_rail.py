"""Tests del rail de herramientas unificado y los flyouts (Fase 4).

OpenSpec: ``2026-09-25-modernize-ui-tool-rail-flyouts``.

Cubre (contra la ventana real, offscreen):
- Paridad 1:1: los 15 botones de categoría del rail y las celdas de flyout
  equivalentes a los ``QMenu``/paletas históricas (2/11/11/10/16/10/12/2/8/23)
  + pies de flyout (anillo/átomo: acción de texto; energía: 2 submenús de
  preset). Clean2D/Validar/Numerar no son botones permanentes del rail
  (convergencia visual con el spike): sus ``QAction`` permanecen accesibles
  por menús y atajos (Ctrl+K).
- Delegación sin lógica propia: el clic en botones/celdas dispara los
  ``QAction``/callbacks originales (canvas + señales de los toolbars); el rail
  no crea ``QAction`` propios;
- Estado activo derivado (``tool_changed`` → highlight; ``clear`` → sin
  highlight; sin emisiones duplicadas de ``tool_changed``);
- Atajos de letra simple contextuales (11 teclas, sin ``QShortcut``; se
  suprimen con modificadores, foco en entrada de texto);
- Toolbars históricas ocultas (no eliminadas) y propietarias de las acciones;
- Métricas (rail 58 px, flyout 244 px) y QSS por tokens (main vs. paletas).
"""

from __future__ import annotations

import pytest
from PyQt6.QtCore import QCoreApplication, QEvent, QPointF, Qt
from PyQt6.QtGui import QKeyEvent, QMouseEvent, QShortcut
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import (
    QApplication,
    QLineEdit,
    QToolButton,
    QWidget,
)

from chemuson.gui.energy_diagrams import (
    ENERGY_DIAGRAM_MENU_ORDER,
    energy_diagram_display_name,
    energy_diagram_tool_id,
)
from chemuson.gui.items import EnergyDiagramItem, TextAnnotationItem
from chemuson.gui.orbitals import (
    ORBITAL_MENU_ORDER,
    orbital_display_name,
    orbital_tool_id,
)
from chemuson.gui.theme import METRICS
from chemuson.gui.theme.qss import (
    get_main_stylesheet,
    get_tool_palette_stylesheet,
)
from chemuson.gui.tool_rail import (
    RAIL_WIDTH,
    ToolRail,
    ToolRailButton,
    ToolShortcutDispatcher,
)
from chemuson.gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar

# ---------------------------------------------------------------------------
# Inventarios de paridad (baseline.md del OpenSpec)
# ---------------------------------------------------------------------------

EXPECTED_RAIL_KEYS = (
    "select",
    "lasso",
    "bond",
    "chain",
    "ring",
    "atom",
    "coord",
    "rotate3d",
    "text",
    "arrows",
    "brackets",
    "symbols",
    "plates",
    "energy",
    "orbitals",
)

EXPECTED_FLYOUT_CELL_COUNTS = {
    "select": 2,  # pointer + lazo
    "bond": 11,  # 3 + 8 tipos de enlace
    "ring": 11,  # benceno + anillos 3..12
    "atom": 10,  # C, N, O, S, P, F, Cl, Br, I, H
    "arrows": 16,  # 16 flechas de anotación
    "brackets": 10,  # corchetes
    "symbols": 12,  # cargas/pares/parcial
    "plates": 2,  # TLC + electroforesis
    "energy": len(ENERGY_DIAGRAM_MENU_ORDER),  # 8 preset
    "orbitals": len(ORBITAL_MENU_ORDER),  # 23 orbitales
}

EXPECTED_SHORTCUT_KEYS = (
    "V",
    "A",
    "L",
    "B",
    "R",
    "C",
    "T",
    "N",
    "G",
    "E",
    "O",
)


def _make_unit_rail():
    """Rail unitario (toolbars reales, sin ventana completa)."""
    window = QWidget()
    window.show()
    toolbar = ChemusonToolbar()
    symbols = SymbolPaletteToolbar(action_group=toolbar.action_group)
    rail = ToolRail(toolbar, symbols, window=window)
    window.layout()
    app = QApplication.instance()
    if app is not None:
        app.processEvents()
    return window, toolbar, symbols, rail


def _foot_button_texts(flyout) -> list[str]:
    texts: list[str] = []
    for i in range(flyout._foot_row.count()):
        item = flyout._foot_row.itemAt(i)
        if item is not None and isinstance(item.widget(), QToolButton):
            texts.append(item.widget().text())
    return texts


# ---------------------------------------------------------------------------
# 1. Paridad 1:1 (botones + celdas + pies)
# ---------------------------------------------------------------------------


def test_rail_button_inventory_is_1_to_1():
    _, _, _, rail = _make_unit_rail()
    assert rail.button_count() == len(EXPECTED_RAIL_KEYS)
    assert sorted(rail._buttons) == sorted(EXPECTED_RAIL_KEYS)
    assert RAIL_WIDTH == METRICS["railW"] == 58


def test_flyout_cell_inventory_matches_historical_menus():
    _, _, _, rail = _make_unit_rail()
    for key, expected in EXPECTED_FLYOUT_CELL_COUNTS.items():
        flyout = rail.open_flyout(key)
        assert flyout is not None, key
        assert len(flyout._cells) == expected, (key, len(flyout._cells))
        flyout.close_with(None)


def test_text_flyout_footer_exposes_color_submenu_in_full_window():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    # El flyout de texto expone el submenú de color como botón de pie
    # (el menú de texto se fija en el ensamblaje con las QAction de la ventana).
    text = rail.open_flyout("text")
    assert _foot_button_texts(text) == ["Colores ▾"]
    text.close_with(None)
    win.close()


def test_flyout_footers_delegate_to_original_actions_and_submenus():
    _, _, _, rail = _make_unit_rail()
    ring = rail.open_flyout("ring")
    assert _foot_button_texts(ring) == ["Tamaño personalizado…"]
    ring.close_with(None)
    atom = rail.open_flyout("atom")
    assert _foot_button_texts(atom) == ["Tabla periódica…"]
    atom.close_with(None)
    energy = rail.open_flyout("energy")
    assert _foot_button_texts(energy) == ["Diagrams ▾", "Presets ▾"]
    energy.close_with(None)


# ---------------------------------------------------------------------------
# 2. Delegación sin lógica propia
# ---------------------------------------------------------------------------


def test_rail_creates_no_own_qactions():
    _, _, _, rail = _make_unit_rail()
    from PyQt6.QtGui import QAction

    assert rail.findChildren(QAction) == []


def test_rail_buttons_delegate_to_original_actions():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    # Cadena: botón del rail -> action_chain del toolbar original.
    rail._buttons["chain"].clicked.emit()
    QApplication.processEvents()
    assert win.canvas.state.active_tool == "tool_chain"
    # Flechas: botón del rail -> action del toolbar de símbolos
    # (emite la flecha actual, por defecto la directa).
    rail._buttons["arrows"].clicked.emit()
    QApplication.processEvents()
    assert win.canvas.state.active_tool == "tool_arrow_forward"
    win.close()


def test_flyout_cell_click_runs_original_toolbar_callback():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    # La celda "Doble" del flyout de enlaces ejecuta el callback del botón del
    # menú original (señal bond_palette_changed del toolbar).
    saw_bond_palette: list[str] = []
    win.toolbar.bond_palette_changed.connect(
        lambda meta: saw_bond_palette.append(meta.get("kind", ""))
    )
    flyout = rail.open_flyout("bond")
    assert flyout is not None
    labels = [cell._label.text().replace("\n", " ").lower() for cell in flyout._cells]
    assert any("doble" in text for text in labels)
    # El orden de celdas replica el menú original (la celda «Enlace doble»
    # es la tercera de la primera fila, con columnas=3).
    index = next(i for i, text in enumerate(labels) if "doble" in text)
    flyout._cells[index].clicked.emit()
    QApplication.processEvents()
    assert saw_bond_palette, "la celda del flyout no ejecutó el callback original"
    win.close()


# ---------------------------------------------------------------------------
# 3. Estado activo derivado (sin duplicar señales)
# ---------------------------------------------------------------------------


def test_active_state_is_derived_from_toolbar_signals():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    rail._buttons["ring"].set_active(False)
    win.toolbar.tool_changed.emit("tool_ring")
    QApplication.processEvents()
    assert rail._buttons["ring"].property("active") is True
    assert rail._buttons["bond"].property("active") is False
    # Normalizado por el canvas (señal de símbolos, tool_id completo).
    win.symbols_toolbar.tool_changed.emit("tool_orbital_s_shaded")
    QApplication.processEvents()
    assert rail._buttons["orbitals"].property("active") is True
    assert rail._buttons["ring"].property("active") is False
    # Cambio de pestaña limpia el highlight.
    win._clear_active_tool_selection()
    QApplication.processEvents()
    assert rail._buttons["orbitals"].property("active") is False
    win.close()


def test_rail_clicks_emit_no_duplicate_tool_changed():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    emissions: list[str] = []
    win.toolbar.tool_changed.connect(emissions.append)
    win.toolbar.tool_changed.emit("tool_select")
    QApplication.processEvents()
    before = len(emissions)
    rail._buttons["bond"].clicked.emit()
    QApplication.processEvents()
    # Un solo nuevo ``tool_changed`` por el disparo del QAction original.
    assert len(emissions) == before + 1
    assert emissions[-1] == "tool_bond"
    win.close()


# ---------------------------------------------------------------------------
# 4. Atajos de letra simple contextuales (sin QShortcut)
# ---------------------------------------------------------------------------


def test_shortcut_map_covers_the_eleven_keys():
    _, _, _, rail = _make_unit_rail()
    mapping = rail.shortcut_map()
    assert len(mapping) == len(EXPECTED_SHORTCUT_KEYS)
    for letter in EXPECTED_SHORTCUT_KEYS:
        assert int(getattr(Qt.Key, f"Key_{letter}")) in mapping


def test_window_registers_no_qshortcut():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.findChildren(QShortcut) == []
    win.close()


def _press_key(target, key: int, modifiers=Qt.KeyboardModifier.NoModifier) -> None:
    event = QKeyEvent(QKeyEvent.Type.KeyPress, key, modifiers)
    QCoreApplication.sendEvent(target, event)
    QApplication.processEvents()


def _reset_tool(win) -> None:
    win.toolbar.tool_changed.emit("tool_select")
    QApplication.processEvents()


def test_shortcuts_activate_tools_and_are_contextual():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()

    def expect(letter: str, expected_tool: str) -> None:
        _reset_tool(win)
        _press_key(win, int(getattr(Qt.Key, f"Key_{letter}")))
        assert win.canvas.state.active_tool == expected_tool

    expect("B", "tool_bond")
    expect("R", "tool_ring")
    expect("L", "tool_chain")
    expect("C", "tool_atom")
    expect("T", "tool_text")
    expect("N", "tool_arrow_forward")
    expect("G", "tool_brackets")
    expect("E", "tool_energy_diagram")
    expect("O", "tool_orbital")
    expect("V", "tool_select")
    expect("A", "tool_select_lasso")
    win.close()


def test_shortcut_with_modifier_is_suppressed():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    _reset_tool(win)
    _press_key(win, int(Qt.Key.Key_B), Qt.KeyboardModifier.ControlModifier)
    assert win.canvas.state.active_tool == "tool_select"
    win.close()


def test_shortcut_suppressed_with_text_input_focus():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    _reset_tool(win)
    line = QLineEdit(win)
    line.show()
    line.setFocus()
    QApplication.processEvents()
    _press_key(line, int(Qt.Key.Key_R))
    assert win.canvas.state.active_tool == "tool_select"
    win.close()


def test_window_builds_without_exception():
    """ChemusonWindow debe poder construirse en todo momento (regresión).

    La ausencia de metadata perfecta en una celda de paleta NO debe ser
    una excepción fatal: los flyouts se reconstruyen con callback estable o
    con el callback histórico como fallback.
    """
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.close()


def test_dispatcher_suppresses_keys_during_graphical_text_and_energy_edit():
    """Las letras de atajo NO deben cambiar de herramienta mientras se
    edita contenido del canvas; el evento sigue su curso normal hacia el
    receptor con foco (el texto llega al item; la edición sigue activa).
    """
    from chemuson.gui.main_window import ChemusonWindow

    # --- Texto: las letras llegan al item sin cambiar herramienta/rail ---
    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    win.toolbar.tool_changed.emit("tool_text")
    item = TextAnnotationItem("", 100.0, 100.0)
    win.canvas.add_text_item(item)
    item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
    item.setFocus()
    win.canvas.remember_text_edit_item(item)
    QApplication.processEvents()
    expected = ""
    for character in ("V", "A", "L"):
        QTest.keyClick(win.canvas.viewport(), character)
        QApplication.processEvents()
        expected += character
        assert item.toPlainText() == expected
        assert item.textInteractionFlags() == Qt.TextInteractionFlag.TextEditorInteraction
        assert win.canvas.current_tool == "tool_text"
        assert win.tool_rail._current_group == "text"
    # Al terminar la edición los atajos vuelven a funcionar.
    win.toolbar.tool_changed.emit("tool_select")
    QApplication.processEvents()
    assert item.textInteractionFlags() == Qt.TextInteractionFlag.NoTextInteraction
    _press_key(win, int(Qt.Key.Key_B))
    assert win.canvas.current_tool == "tool_bond"
    win.canvas.undo_stack.clear()
    win.canvas.undo_stack.setClean()
    win.close()

    # --- Energía: edición directa no roba las letras ni cambia herramienta ---
    win2 = ChemusonWindow()
    win2.show()
    QApplication.processEvents()
    win2.canvas.state.active_energy_diagram_kind = "sublevel_p"
    energy_item = win2.canvas._insert_energy_diagram_item(QPointF(120.0, 120.0))
    energy_item.begin_direct_edit()
    QApplication.processEvents()
    assert isinstance(energy_item, EnergyDiagramItem)
    assert energy_item.is_editing()
    _press_key(win2, int(Qt.Key.Key_C))
    assert win2.canvas.current_tool != "tool_atom"
    assert energy_item.is_editing()
    energy_item.end_direct_edit()
    _press_key(win2, int(Qt.Key.Key_C))
    assert win2.canvas.current_tool == "tool_atom"
    win2.canvas.undo_stack.clear()
    win2.canvas.undo_stack.setClean()
    win2.close()


def test_dispatcher_suppresses_keys_outside_the_window():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    foreign = QWidget()
    foreign.show()
    _reset_tool(win)
    _press_key(foreign, int(Qt.Key.Key_R))
    assert win.canvas.state.active_tool == "tool_select"
    win.close()
    foreign.close()


# ---------------------------------------------------------------------------
# 5. Toolbars históricas: ocultas, no eliminadas, propietarias
# ---------------------------------------------------------------------------


def test_historical_toolbars_are_hidden_but_still_own_actions():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.toolbar.isVisible() is False
    assert win.symbols_toolbar.isVisible() is False
    # Las acciones siguen vivas y propietarias del grupo exclusivo.
    assert win.toolbar.action_group.actions()
    assert win.toolbar.select_action in win.toolbar.action_group.actions()
    assert win.canvas.set_current_tool is not None  # handlers intactos
    win.close()


# ---------------------------------------------------------------------------
# 6. Métricas + QSS por tokens
# ---------------------------------------------------------------------------


def test_rail_and_flyout_metrics_match_mockup():
    assert METRICS["railW"] == 58
    assert METRICS["flyoutW"] == 244
    # Métricas del spike aprobado (tokens espejo del spike).
    assert METRICS["railBtn"] == 42
    assert METRICS["railIcon"] == 21
    assert METRICS["statusH"] == 34
    assert METRICS["appbarH"] == 54
    _, _, _, rail = _make_unit_rail()
    assert rail.width() == 58
    flyout = rail.open_flyout("bond")
    assert flyout.width() == 244
    flyout.close_with(None)


# ---------------------------------------------------------------------------
# 6b. Convergencia visual con el spike aprobado
# ---------------------------------------------------------------------------


def test_shell_has_no_visible_menubar_or_toolbars():
    """La QMenuBar y las toolbars clásicas no se muestran; la QMenuBar
    sigue existiendo (menús/acciones/atajos intactos) y se abre como popup
    (hamburguesa o Alt)."""
    from PyQt6.QtWidgets import QToolBar

    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.menuBar() is not None
    assert win.menuBar().isVisible() is False
    # Los 6 menús superiores siguen en el menubar oculto.
    assert len([a for a in win.menuBar().actions() if a.menu() is not None]) == 6
    visible = [t.objectName() for t in win.findChildren(QToolBar) if t.isVisible()]
    assert visible == []
    win.close()


def test_rail_is_bare_widget_in_central_layout():
    """El rail es un ``QWidget`` de 58 px dentro del layout central (no un
    contenedor ``QToolBar``)."""
    from PyQt6.QtWidgets import QToolBar

    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    assert not isinstance(rail, QToolBar)
    # No está emparentado a ningún QToolBar.
    parent = rail.parentWidget()
    assert parent is None or not isinstance(parent, QToolBar)
    win.close()


def test_second_click_on_active_category_opens_flyout():
    """Clic izquierdo: primer clic activa; un segundo clic sobre la misma
    categoría activa abre su flyout."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    win.toolbar.tool_changed.emit("tool_select")
    QApplication.processEvents()
    # Primer clic: activa.
    rail._buttons["bond"].clicked.emit()
    QApplication.processEvents()
    assert win.canvas.state.active_tool == "tool_bond"
    # Segundo clic sobre la categoría activa: flyout.
    rail._buttons["bond"].clicked.emit()
    QApplication.processEvents()
    assert rail._flyouts["bond"].isVisible()
    rail._flyouts["bond"].close_with(None)
    # Una categoría inactiva con menú no abre flyout al activarse.
    rail._buttons["atom"].clicked.emit()
    QApplication.processEvents()
    assert win.canvas.state.active_tool == "tool_atom"
    assert not rail._flyouts["atom"].isVisible()
    win.close()


def test_text_toolbar_hidden_by_default_and_contextual():
    """La toolbar de texto está oculta por defecto y aparece solo con la
    herramienta de texto (sus acciones/señales no cambian)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.text_toolbar.isVisible() is False
    win.toolbar.tool_changed.emit("tool_text")
    QApplication.processEvents()
    assert win.text_toolbar.isVisible() is True
    win.toolbar.tool_changed.emit("tool_select")
    QApplication.processEvents()
    assert win.text_toolbar.isVisible() is False
    win.close()


def test_shell_metrics_match_spike():
    """Alturas del shell: appbar 54 px, rail 58 px, status 34 px; la
    hamburguesa existe con icono y tamaño mínimo 900×560."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.resize(1440, 900)
    win.show()
    QApplication.processEvents()
    assert win.app_bar.height() == 54
    assert win.tool_rail.width() == 58
    assert win.statusBar().height() == 34
    assert win.app_bar.menu_button is not None
    assert not win.app_bar.menu_button.icon().isNull()
    assert (win.minimumWidth(), win.minimumHeight()) == (900, 560)
    win.close()


def test_window_fits_980x600():
    """El shell completo cabe en 980×600 (el rail se desplaza compacto
    con scroll invisible; sin micro-iconos)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.resize(980, 600)
    win.show()
    QApplication.processEvents()
    assert win.width() <= 980
    assert win.height() <= 600
    # Los botones conservan el tamaño del spike (42 px) aunque el rail
    # precise scroll.
    assert win.tool_rail._buttons["select"].width() == 42
    win.close()


def test_clean2d_validate_numbering_remain_accessible_without_rail_buttons():
    """Las tres acciones salen del rail permanente pero siguen vivas:
    ``QAction`` existentes, menús que las contienen y atajo Ctrl+K."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    def in_some_menu(action) -> bool:
        for top in win.menuBar().actions():
            top_menu = top.menu()
            if top_menu is None:
                continue
            for sub in top_menu.actions():
                holder = sub.menu() if sub.menu() is not None else top_menu
                if action in holder.actions():
                    return True
        return False

    for attr in ("action_clean_2d_full", "action_validate_structure",
                 "action_numbering_recalculate"):
        action = getattr(win, attr)
        assert action is not None
        # Cada una vive en algún QMenu (menús intactos del menubar oculto).
        assert in_some_menu(action), attr
    assert win.action_clean_2d_full.shortcut().toString() == "Ctrl+K"
    win.close()


def test_alt_opens_menu_popup_and_hamburger_exists():
    """Alt (sin modificadores) y la hamburgensa de la app bar no rompen el
    shell (el popup de la QMenuBar se ancla bajo la hamburgensa)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    _press_key(win, int(Qt.Key.Key_Alt))
    win.app_bar.menu_button.clicked.emit()
    QApplication.processEvents()
    # El menubar oculto sigue siendo la fuente de los 6 menús.
    assert len([a for a in win.menuBar().actions() if a.menu() is not None]) == 6
    win.close()


def test_no_kbd_badges_in_rail():
    """Los badges kbd visibles se eliminan del rail normal; el atajo se
    documenta en el tooltip y sigue funcional."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    for key, button in win.tool_rail._buttons.items():
        assert button._kbd is None, key
    # El tooltip conserva la pista de atajo (p. ej. "Enlaces (B)").
    assert "(B)" in win.tool_rail._buttons["bond"].toolTip()
    win.close()


def test_rail_and_flyout_qss_use_tokens_in_main_stylesheet_only():
    for theme in ("light", "dark"):
        main_qss = get_main_stylesheet(theme)
        for selector in (
            "#toolRail",
            "QToolButton#railBtn",
            'QToolButton#railBtn[active="true"]',
            "#railKbd",
            "#railSep",
            "#flyout",
            "#flyoutTitle",
            '#flyKbd',
            'QFrame[cls="flyItem"]',
            'QToolButton[cls="flyFoot"]',
            "#flyFootTxt",
        ):
            assert selector in main_qss, (theme, selector)
        # El QSS de paletas históricas no recibe estilos del rail.
        palette_qss = get_tool_palette_stylesheet(theme)
        assert "#toolRail" not in palette_qss
        assert "#flyout" not in palette_qss


# ---------------------------------------------------------------------------
# 7. Refresh de tema (iconos + flyouts reconstruidos)
# ---------------------------------------------------------------------------


def test_theme_refresh_rebuilds_rail_and_flyouts():
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    rail = win.tool_rail
    counts_before = {
        key: len(rail.open_flyout(key)._cells)
        for key in ("bond", "ring", "energy", "orbitals")
    }
    for flyout in list(rail._flyouts.values()):
        flyout.close_with(None)
    win.current_theme = "dark"
    win._apply_theme()
    QApplication.processEvents()
    for key, expected in counts_before.items():
        flyout = rail.open_flyout(key)
        assert len(flyout._cells) == expected, key
        flyout.close_with(None)
    win.current_theme = "light"
    win._apply_theme()
    QApplication.processEvents()
    win.close()


def test_button_widget_api():
    button = ToolRailButton()
    assert button.objectName() == "railBtn"
    button.set_active(True)
    assert button.property("active") is True
    button.set_kbd("B")
    assert button._kbd is not None
    button.set_kbd(None)
    assert button._kbd is None


@pytest.mark.parametrize(
    ("tool_id", "expected"),
    [
        ("tool_select", "select"),
        ("tool_select_lasso", "select"),
        ("bond_double", "bond"),
        ("tool_ring", "ring"),
        ("atom_c", "atom"),
        ("tool_chain", "chain"),
        ("coord_5", "coord"),
        ("tool_rotate_3d_precise", "rotate3d"),
        ("tool_text", "text"),
        ("tool_arrow_up", "arrows"),
        ("tool_brackets_round", "brackets"),
        ("tool_tlc", "plates"),
        ("tool_charge_plus", "symbols"),
        ("tool_energy_diagram_sublevel_s", "energy"),
        ("tool_orbital_s_shaded", "orbitals"),
        ("tool_none", None),
        ("", None),
    ],
)
def test_group_for_tool_normalization(tool_id, expected):
    from chemuson.gui.tool_rail import _group_for_tool

    assert _group_for_tool(tool_id) == expected


def test_shortcut_dispatcher_class_exists():
    assert ToolShortcutDispatcher is not None


def test_ctrl_k_still_triggers_clean_2d_full_and_menubar_hidden():
    """El dispatcher no roba combinaciones con modificadores: Ctrl+K sigue
    perteneciendo a ``action_clean_2d_full`` (Fase 3). La QMenuBar queda
    oculta (convergencia con el spike) pero viva: sus QMenu/QAction/
    atajos son los originales y se abren por popup (hamburguesa/Alt)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.menuBar() is not None
    assert win.menuBar().isVisible() is False

    saw_clean_2d = []
    win.action_clean_2d_full.triggered.connect(lambda: saw_clean_2d.append(True))
    _reset_tool(win)
    _press_key(win, int(Qt.Key.Key_K), Qt.KeyboardModifier.ControlModifier)
    # Ctrl+K sigue disparando la acción original (el dispatcher lo ignora:
    # modificadores presentes).
    assert saw_clean_2d == [True], "Ctrl+K debe seguir disparando action_clean_2d_full"
    # Y no activó la herramienta K (no existe) ni cambió el grupo activo.
    assert win.canvas.state.active_tool == "tool_select"
    assert win.tool_rail._current_group == "select"
    win.close()


# ---------------------------------------------------------------------------
# Reconexion de las paletas avanzadas del rail (orbitales, energía,
# símbolos, corchetes) + regresión del crash de flyout (globalPos)
# ---------------------------------------------------------------------------


def _flyout_cell_by_tooltip(flyout, tooltip: str):
    for cell in flyout._cells:
        if cell.toolTip() == tooltip:
            return cell
    return None


def test_flyout_outside_click_does_not_raise_and_closes():
    """Regresión de crash_20260926_080241/080409: PyQt6 removió
    ``QMouseEvent.globalPos()``; el filtro de cierre por clic fuera debe
    usar ``globalPosition().toPoint()`` y cerrar el flyout sin lanzar."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    flyout = win.tool_rail.open_flyout("orbitals")
    assert flyout is not None
    assert flyout.isVisible()
    # Clic fuera (esquina inferior derecha, siempre fuera del flyout por el
    # clamp de 8 px) enviado a la ventana: el filtro debe ejecutar sin
    # AttributeError y cerrar el flyout.
    ev = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(win.width() - 2.0, win.height() - 2.0),
        QPointF(win.width() - 2.0, win.height() - 2.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    QApplication.sendEvent(win, ev)
    QApplication.processEvents()
    assert flyout.isVisible() is False
    win.canvas.undo_stack.clear()
    win.canvas.undo_stack.setClean()
    win.close()


def test_orbital_end_to_end_from_rail_flyout():
    """Rail → flyout → clic en celda → herramienta orbital activa → clic en
    canvas crea el orbital, con undo/redo."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        flyout = win.tool_rail.open_flyout("orbitals")
        assert flyout is not None
        assert flyout.isVisible()
        kind = "sigma_bonding_solid"
        cell = _flyout_cell_by_tooltip(flyout, orbital_display_name(kind))
        assert cell is not None, "El flyout debe exponer la celda sigma_bonding_solid"
        QTest.mouseClick(cell, Qt.MouseButton.LeftButton)
        QApplication.processEvents()
        assert win._current_tool_id == orbital_tool_id(kind)
        assert win.canvas.current_tool == "tool_orbital"
        assert win.canvas.state.active_orbital_kind == kind
        assert win.symbols_toolbar.orbital_action.isChecked() is True
        # El flyout se cierra tras la selección.
        assert flyout.isVisible() is False
        # El clic de canvas crea el item con el kind elegido.
        before = len(win.canvas.orbital_items)
        win.canvas._insert_orbital_item(QPointF(150.0, 150.0))
        assert len(win.canvas.orbital_items) == before + 1
        assert win.canvas.orbital_items[-1].kind() == kind
        win.canvas.undo_stack.undo()
        assert len(win.canvas.orbital_items) == before
        win.canvas.undo_stack.redo()
        assert len(win.canvas.orbital_items) == before + 1
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def test_energy_diagram_end_to_end_from_rail_flyout():
    """Rail → flyout → clic en celda → herramienta de diagrama activa → clic
    en canvas crea el diagrama, con undo/redo."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        flyout = win.tool_rail.open_flyout("energy")
        assert flyout is not None
        assert flyout.isVisible()
        kind = "hybrid_sp2"
        cell = _flyout_cell_by_tooltip(flyout, energy_diagram_display_name(kind))
        assert cell is not None, "El flyout debe exponer la celda hybrid_sp2"
        QTest.mouseClick(cell, Qt.MouseButton.LeftButton)
        QApplication.processEvents()
        assert win._current_tool_id == energy_diagram_tool_id(kind)
        assert win.canvas.current_tool == "tool_energy_diagram"
        assert win.canvas.state.active_energy_diagram_kind == kind
        assert win.symbols_toolbar.energy_diagram_action.isChecked() is True
        assert flyout.isVisible() is False
        before = len(
            [item for item in win.canvas.scene.items() if isinstance(item, EnergyDiagramItem)]
        )
        win.canvas._insert_energy_diagram_item(QPointF(150.0, 150.0))
        after = len(
            [item for item in win.canvas.scene.items() if isinstance(item, EnergyDiagramItem)]
        )
        assert after == before + 1
        win.canvas.undo_stack.undo()
        assert (
            len([item for item in win.canvas.scene.items() if isinstance(item, EnergyDiagramItem)])
            == before
        )
        win.canvas.undo_stack.redo()
        assert (
            len([item for item in win.canvas.scene.items() if isinstance(item, EnergyDiagramItem)])
            == after
        )
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def test_symbol_cell_click_selects_symbol_tool():
    """El clic en una celda de símbolos del flyout selecciona el símbolo
    original (delegación 1:1, sin crash)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        flyout = win.tool_rail.open_flyout("symbols")
        assert flyout is not None
        cell = _flyout_cell_by_tooltip(flyout, "Carga positiva")
        assert cell is not None
        QTest.mouseClick(cell, Qt.MouseButton.LeftButton)
        QApplication.processEvents()
        assert win._current_tool_id == "tool_charge_plus"
        assert win.canvas.current_tool == "tool_charge_plus"
        assert win.symbols_toolbar.symbol_action.isChecked() is True
        assert win.tool_rail._current_group == "symbols"
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def test_bracket_cell_click_selects_bracket_tool():
    """Los 10 corchetes históricos están en el flyout y el clic selecciona
    la herramienta de corchetes original (sin RuntimeError de callback)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        flyout = win.tool_rail.open_flyout("brackets")
        assert flyout is not None
        assert len(flyout._cells) == 10
        QTest.mouseClick(flyout._cells[0], Qt.MouseButton.LeftButton)  # tool_brackets_square (orden histórico)
        QApplication.processEvents()
        assert win.canvas.current_tool == "tool_brackets"
        assert win.canvas.state.active_bracket_type == "[]"
        assert win.symbols_toolbar.bracket_action.isChecked() is True
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def test_energy_submenus_and_presets_connected():
    """Los submenús de energía (diagramas electrónicos y presets) conservan
    su wiring original: la acción del botón emite ``tool_changed``; las
    QActions de los submenús están conectadas a las señales de solicitud y
    esas señales tienen el handler de la ventana (sin disparar los diálogos
    modales: se verifica la conexión, no la ejecución). El flyout de energía
    expone los dos submenús como botones de pie."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        # Botón (acción) a nivel de toolbar: emite la herramienta actual.
        win.symbols_toolbar.energy_diagram_action.trigger()
        QApplication.processEvents()
        assert win.canvas.current_tool == "tool_energy_diagram"

        # Submenús vivos dentro del menú de energía.
        menu = win.symbols_toolbar.energy_diagram_button.menu()
        submenus = [action.menu() for action in menu.actions() if action.menu() is not None]
        assert len(submenus) == 2, "Faltan los submenús de energía"
        electronic_menu, presets_menu = submenus

        # Wiring a nivel de toolbar: la QAction emite la señal de solicitud.
        atomic_action = electronic_menu.actions()[0]
        assert atomic_action.receivers(atomic_action.triggered) > 0
        toolbar = win.symbols_toolbar
        assert toolbar.receivers(toolbar.atomic_diagram_requested) >= 1
        assert toolbar.receivers(toolbar.diatomic_mo_diagram_requested) >= 1
        assert toolbar.receivers(toolbar.ligand_field_diagram_requested) >= 1

        # Wiring a nivel de ventana: el handler existe (señal con receptores).
        assert win.symbols_toolbar.receivers(win.symbols_toolbar.atomic_diagram_requested) >= 1
        assert (
            win.symbols_toolbar.receivers(
                win.symbols_toolbar.electronic_diagram_preset_requested
            )
            >= 1
        )

        # Presets: tres submenús, cada QAction conectada.
        preset_submenus = [
            action.menu() for action in presets_menu.actions() if action.menu() is not None
        ]
        assert len(preset_submenus) == 3, "Faltan los submenús de presets"
        for preset_menu in preset_submenus:
            assert len(preset_menu.actions()) > 0
            for action in preset_menu.actions():
                assert action.receivers(action.triggered) > 0

        # El flyout de energía expone los dos submenús como botones de pie.
        flyout = win.tool_rail.open_flyout("energy")
        assert flyout is not None
        assert len(flyout._foot_buttons) == 2
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()


def test_theme_refresh_cycle_preserves_flyout_callbacks():
    """light → dark → light no destruye los callbacks del flyout: tras el
    ciclo, el clic en una celda sigue seleccionando la herramienta
    original."""
    from chemuson.gui.main_window import ChemusonWindow
    from chemuson.gui.theme import apply_theme

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        apply_theme(win, "dark")
        QApplication.processEvents()
        apply_theme(win, "light")
        QApplication.processEvents()

        flyout = win.tool_rail.open_flyout("orbitals")
        assert flyout is not None
        cell = _flyout_cell_by_tooltip(flyout, orbital_display_name("dz2_shaded"))
        assert cell is not None
        QTest.mouseClick(cell, Qt.MouseButton.LeftButton)
        QApplication.processEvents()
        assert win.canvas.current_tool == "tool_orbital"
        assert win.canvas.state.active_orbital_kind == "dz2_shaded"
    finally:
        win.canvas.undo_stack.clear()
        win.canvas.undo_stack.setClean()
        win.close()
