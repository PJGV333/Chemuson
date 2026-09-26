"""Tests del rail de herramientas unificado y los flyouts (Fase 4).

OpenSpec: ``2026-09-25-modernize-ui-tool-rail-flyouts``.

Cubre (contra la ventana real, offscreen):
- Paridad 1:1: los 18 botones del rail y las celdas de flyout equivalentes a
  los ``QMenu``/paletas históricas (2/11/11/10/16/10/12/2/8/23) + pies de
  flyout (anillo/átomo: acción de texto; energía: 2 submenús de preset);
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
from PyQt6.QtCore import QCoreApplication, Qt
from PyQt6.QtGui import QKeyEvent, QShortcut
from PyQt6.QtWidgets import (
    QApplication,
    QLineEdit,
    QToolButton,
    QWidget,
)

from chemuson.gui.energy_diagrams import ENERGY_DIAGRAM_MENU_ORDER
from chemuson.gui.orbitals import ORBITAL_MENU_ORDER
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
    "clean2d",
    "validate",
    "numbering",
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
    _, _, _, rail = _make_unit_rail()
    assert rail.width() == 58
    flyout = rail.open_flyout("bond")
    assert flyout.width() == 244
    flyout.close_with(None)


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


def test_ctrl_k_still_triggers_clean_2d_full_and_menubar_visible():
    """El dispatcher no roba combinaciones con modificadores: Ctrl+K sigue
    perteneciendo a ``action_clean_2d_full`` (Fase 3) y la QMenuBar se
    mantiene visible (el rail no la reemplaza)."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    assert win.menuBar().isVisible()

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
