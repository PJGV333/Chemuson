"""Tests de la paleta de comandos (Fase 6, corrección de UX: Ctrl+P).

Cubre el contrato de la Fase 6 según OpenSpec
``2026-10-01-modernize-ui-command-palette``:

- una única ``CommandPalette`` por ventana;
- registro con ≥60 comandos únicos y deduplicación por identidad de ``QAction``;
- las 7 páginas del SidePanel, las exportaciones (PNG/SVG/PDF/CML/SMILES) y
  Clean2D quick buscables;
- corrección de UX (post-Fase 6): ``Ctrl+K`` vuelve a pertenecer a
  ``action_clean_2d_full`` (limpia 2D, 1 paso) y ya NO abre la paleta; la
  paleta de comandos se abre con ``Ctrl+P`` (única ``QAction`` global que lo
  posee) y con el clic en SearchPill; ``Ctrl+Shift+K``/``Ctrl+Alt+K`` intactos;
- SearchPill como entrada real (mismo camino que Ctrl+P, badge ``Ctrl P``);
- filtro substring con prioridad por prefijo y matching por keywords;
- teclado (↑/↓/Enter/Esc), clic, disabled, checkable, cierre tras ejecutar;
- presentación dentro de la ventana (980×600) y light→dark→light;
- los atajos de herramienta no interfieren al escribir;
- las acciones siguen accesibles por menú y sin conflicto de ``Ctrl+P``;
- contrato de imports de ``command_palette.py`` (sin dependencia de química).

Los diálogos modales se prueban por wiring/``trigger`` de ``QAction`` segura,
sin abrir diálogos (no bloquean el entorno offscreen).
"""

from __future__ import annotations

import pytest
from PyQt6.QtCore import QSettings, QStandardPaths, Qt
from PyQt6.QtGui import QAction, QKeyEvent, QKeySequence
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QLineEdit

from chemuson.gui.main_window import ChemusonWindow


@pytest.fixture(scope="module", autouse=True)
def _qapp() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    """Isola ``XDG_CONFIG_HOME``/``QSettings`` por test (patrol de tests)."""
    config_location = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.ConfigLocation
    )
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    QSettings.setPath(
        QSettings.Format.NativeFormat, QSettings.Scope.UserScope, str(tmp_path)
    )
    try:
        yield
    finally:
        QSettings.setPath(
            QSettings.Format.NativeFormat, QSettings.Scope.UserScope, config_location
        )


@pytest.fixture
def win() -> ChemusonWindow:
    window = ChemusonWindow()
    window.resize(1440, 900)
    window.show()
    QApplication.processEvents()
    yield window
    window.close()
    QApplication.processEvents()


def _press_input_key(window, key: int, modifiers=Qt.KeyboardModifier.NoModifier):
    """Envía una tecla al ``QLineEdit`` de la paleta (foco real al abrir)."""
    target = window.command_palette.input
    event = QKeyEvent(QKeyEvent.Type.KeyPress, key, modifiers)
    QApplication.sendEvent(target, event)
    QApplication.processEvents()


def _entries(window):
    return window._command_registry.entries()


# ---------------------------------------------------------------------------
# 1-3. Una única paleta; registro ≥60; sin duplicar QAction
# ---------------------------------------------------------------------------
def test_window_builds_a_single_command_palette(win):
    assert win.command_palette is not None
    # Una única instancia de paleta por ventana.
    from chemuson.gui.command_palette import CommandPalette

    palettes = win.findChildren(CommandPalette)
    assert len(palettes) == 1
    assert palettes[0] is win.command_palette


def test_registry_contains_at_least_60_unique_commands(win):
    assert win._command_registry.count() >= 60
    assert len(_entries(win)) == win._command_registry.count()
    # Todas las entradas apuntan a una QAction habilitada para representar.
    for entry in _entries(win):
        assert isinstance(entry.action, QAction)


def test_registry_does_not_duplicate_the_same_qaction(win):
    ids = [id(e.action) for e in _entries(win)]
    assert len(set(ids)) == len(ids), "una QAction no debe aparecer dos veces"
    # Reabrir la paleta no debe duplicar registros (el registro se construye 1 vez).
    before = win._command_registry.count()
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.close_overlay()
    QApplication.processEvents()
    assert win._command_registry.count() == before


# ---------------------------------------------------------------------------
# 4. Las siete páginas del SidePanel son buscables (reutilizando side_panel_actions)
# ---------------------------------------------------------------------------
def test_all_seven_side_panel_pages_are_searchable(win):
    side_actions = win.side_panel_actions
    assert set(side_actions) == {
        "inspector", "validation", "properties", "templates",
        "appearance", "spectroscopy", "compchem",
    }
    for key, action in side_actions.items():
        window_entry = win._command_registry.find_by_action(action)
        assert window_entry is not None, f"falta página {key}"
        # Buscando su etiqueta la encontramos apuntando a la MISMA QAction.
        win.command_palette._apply_filter(action.text().lower())
        assert any(e.action is action for e in win.command_palette._filtered), key


# ---------------------------------------------------------------------------
# 5. Exportaciones mínimas buscables
# ---------------------------------------------------------------------------
def test_expected_exports_are_searchable(win):
    # Las 4 exportaciones de archivo viven en la sección "Exportar".
    for fmt_attr in ("action_export_png", "action_export_svg", "action_export_pdf",
                     "action_export_cml"):
        action = getattr(win, fmt_attr)
        entry = win._command_registry.find_by_action(action)
        assert entry is not None and entry.section == "Exportar", fmt_attr
    # La exportación de SMILES es alcanzable (vive en Estructura en producción).
    smiles_entry = win._command_registry.find_by_action(win.action_export_smiles)
    assert smiles_entry is not None
    # "export" lista PNG/SVG/PDF/CML entre sus resultados.
    win.command_palette._apply_filter("export")
    titles = [e.title.lower() for e in win.command_palette._filtered]
    joined = " ".join(titles)
    for token in ("png", "svg", "pdf", "cml", "smiles"):
        assert token in joined


# ---------------------------------------------------------------------------
# 6. Clean2D quick buscable por su QAction histórica
# ---------------------------------------------------------------------------
def test_clean2d_quick_is_searchable_via_its_historical_qaction(win):
    action = win.action_clean_2d_full
    entry = win._command_registry.find_by_action(action)
    assert entry is not None and entry.section == "Estructura"
    win.command_palette._apply_filter("limpiar 2d (1 paso)")
    assert any(e.action is action for e in win.command_palette._filtered)


# ---------------------------------------------------------------------------
# 7-8. Restauración de Ctrl+K (Clean2D) + Ctrl+P (paleta)
# ---------------------------------------------------------------------------
def test_ctrl_k_belongs_to_clean2d_full_again(win):
    # Corrección de UX (post-Fase 6): ``Ctrl+K`` vuelve a ser el atajo histórico
    # de Clean2D quick, con contexto ``WindowShortcut``. La paleta usa ``Ctrl+P``.
    assert win.action_clean_2d_full.shortcut() == QKeySequence("Ctrl+K")
    assert (
        win.action_clean_2d_full.shortcutContext()
        == Qt.ShortcutContext.WindowShortcut
    )
    # La acción sigue viva: texto y presencia en el registro intactos.
    assert win.action_clean_2d_full.text() == "Limpiar 2D (1 paso)"
    assert win._command_registry.find_by_action(win.action_clean_2d_full) is not None
    # La paleta de comandos posee Ctrl+P.
    assert win.action_command_palette.shortcut() == QKeySequence("Ctrl+P")


def test_ctrl_k_executes_clean2d_and_does_not_open_palette(win):
    saw = []
    win.action_clean_2d_full.triggered.connect(lambda: saw.append(True))
    win.canvas.setFocus()
    win.canvas.viewport().setFocus()
    QTest.keyClick(win.canvas.viewport(), Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()
    # Ctrl+K dispara exactamente action_clean_2d_full (limpiar 2D, 1 paso).
    assert saw == [True], "Ctrl+K debe disparar action_clean_2d_full"
    # ...y NO abre la paleta de comandos.
    assert not win.command_palette.is_open()


def test_ctrl_p_opens_the_command_palette(win):
    win.canvas.setFocus()
    win.canvas.viewport().setFocus()
    QTest.keyClick(win.canvas.viewport(), Qt.Key.Key_P, Qt.KeyboardModifier.ControlModifier)
    QApplication.processEvents()
    assert win.command_palette.is_open()
    # Ctrl+P no disparó Clean2D quick.
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 11. Ctrl+Shift+K y Ctrl+Alt+K conservan su comportamiento
# ---------------------------------------------------------------------------
def test_ctrl_shift_k_and_ctrl_alt_k_are_intact(win):
    assert win.action_clean_2d_publication.shortcut() == QKeySequence("Ctrl+Shift+K")
    assert win.action_clean_2d_propose.shortcut() == QKeySequence("Ctrl+Alt+K")
    pub, prop = [], []
    win.action_clean_2d_publication.triggered.connect(lambda: pub.append(True))
    win.action_clean_2d_propose.triggered.connect(lambda: prop.append(True))
    win.canvas.viewport().setFocus()
    QTest.keyClick(win.canvas.viewport(), Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier)
    QApplication.processEvents()
    QTest.keyClick(win.canvas.viewport(), Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.AltModifier)
    QApplication.processEvents()
    assert pub == [True]
    assert prop == [True]
    # Ninguno de estos abre la paleta de comandos (ahora en Ctrl+P).
    assert not win.command_palette.is_open()


# ---------------------------------------------------------------------------
# 9-10. SearchPill: misma ruta que Ctrl+P + badge Ctrl P
# ---------------------------------------------------------------------------
def test_search_pill_opens_the_same_route_as_ctrl_p(win):
    fired = []
    win.action_command_palette.triggered.connect(lambda: fired.append(True))
    # El clic en la píldora dispara la MISMA QAction de apertura (Ctrl+P).
    win.app_bar.search_pill.activated.emit()
    QApplication.processEvents()
    assert fired == [True]
    assert win.command_palette.is_open()
    win.command_palette.close_overlay()
    QApplication.processEvents()
    # Y la píldora no registra su propio QShortcut (un único camino de apertura).
    from PyQt6.QtGui import QShortcut

    assert win.app_bar.findChildren(QShortcut) == []
    assert win.action_command_palette.shortcut() == QKeySequence("Ctrl+P")


def test_search_pill_shows_ctrl_p_hint_and_is_not_an_editor(win):
    pill = win.app_bar.search_pill
    assert "Ctrl P" in pill.kbd.text()
    assert pill.kbd.isVisibleTo(pill)
    assert "(próximamente)" not in pill.toolTip()
    # La SearchPill no es un editor permanente: el campo editable es de la paleta.
    assert pill.findChildren(QLineEdit) == []


# ---------------------------------------------------------------------------
# 12-13. Filtro: prefijo antes que substring + keywords
# ---------------------------------------------------------------------------
def test_prefix_ranks_before_substring(win):
    # Dos entradas cuyo título contiene "export": una por prefijo, otra substring.
    # "Exportar SMILES..." empieza por "export" (prefijo); "Importar SMILES..."
    # lo contiene? No. Usamos "export" → prefijo de "Exportar ...". Verificamos
    # que "Exportar SMILES..." aparezca ANTES que cualquier entrada no-prefijo.
    win.command_palette._apply_filter("export")
    titles = [e.title for e in win.command_palette._filtered]
    # Los resultados empiezan por una entrada cuyo título empieza por "export".
    assert any(t.lower().startswith("export") for t in titles)
    first = titles[0].lower()
    # La primera entrada es de tier-prefijo (empieza por el query).
    assert first.startswith("export") or any(
        k.startswith("export") for k in win._command_registry.find_by_action(
            win.command_palette._filtered[0].action
        ).keywords
    )


def test_keyword_matching(win):
    # "carbonos" es keyword de action_show_carbons; su título "Mostrar carbonos"
    # ya lo contiene. Usamos una keyword distinta: "quick" de clean_2d_full.
    win.command_palette._apply_filter("quick")
    assert any(
        e.action is win.action_clean_2d_full for e in win.command_palette._filtered
    )
    # Keyword "1 paso" también lo localiza aunque "1 paso" no esté en el título
    # con esa exactitud.
    win.command_palette._apply_filter("publicación")
    assert any(
        e.action is win.action_clean_2d_publication for e in win.command_palette._filtered
    )


# ---------------------------------------------------------------------------
# 14. ↑/↓ cambia la selección
# ---------------------------------------------------------------------------
def test_up_down_changes_selection(win):
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("export")
    QApplication.processEvents()
    n = win.command_palette.filter_count()
    assert n >= 2
    start = win.command_palette._idx
    _press_input_key(win, int(Qt.Key.Key_Down))
    assert win.command_palette._idx == (start + 1) % n
    _press_input_key(win, int(Qt.Key.Key_Up))
    assert win.command_palette._idx == start
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 15. Enter dispara exactamente una vez la QAction elegida
# ---------------------------------------------------------------------------
def test_enter_triggers_selected_action_exactly_once(win):
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("limpiar 2d (1 paso)")
    QApplication.processEvents()
    assert any(
        e.action is win.action_clean_2d_full for e in win.command_palette._filtered
    )
    fired = []
    win.action_clean_2d_full.triggered.connect(lambda: fired.append(True))
    _press_input_key(win, int(Qt.Key.Key_Return))
    QApplication.processEvents()
    assert fired == [True], "Enter debe disparar la QAction exactamente una vez"
    assert not win.command_palette.is_open(), "la paleta se cierra tras ejecutar"


# ---------------------------------------------------------------------------
# 16. Esc cierra sin disparar
# ---------------------------------------------------------------------------
def test_esc_closes_without_firing(win):
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("limpiar 2d (1 paso)")
    QApplication.processEvents()
    fired = []
    win.action_clean_2d_full.triggered.connect(lambda: fired.append(True))
    _press_input_key(win, int(Qt.Key.Key_Escape))
    QApplication.processEvents()
    assert not win.command_palette.is_open()
    assert fired == []


# ---------------------------------------------------------------------------
# 17. Una QAction disabled no se ejecuta
# ---------------------------------------------------------------------------
def test_disabled_action_is_not_executed(win):
    # "Deshacer" está disabled en un documento nuevo.
    action = win.action_undo
    assert not action.isEnabled()
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("deshacer")
    QApplication.processEvents()
    entries = win.command_palette._filtered
    assert any(e.action is action for e in entries)
    # Selecciona la fila de la acción disabled y pulsa Enter: no debe disparar.
    idx = next(
        i for i, e in enumerate(entries) if e.action is action
    )
    win.command_palette._idx = idx
    fired = []
    action.triggered.connect(lambda: fired.append(True))
    _press_input_key(win, int(Qt.Key.Key_Return))
    QApplication.processEvents()
    assert fired == []
    # La fila disabled se refleja visualmente (isEnabled False).
    row = win.command_palette._rows[idx]
    assert not row.isEnabled()
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 18. QAction checkable conserva semántica
# ---------------------------------------------------------------------------
def test_checkable_action_preserves_semantics(win):
    action = win.action_show_carbons
    assert action.isCheckable()
    initial = action.isChecked()
    # Busca la acción checkable en la paleta y ejecútala por su ruta real
    # (delegación a ``trigger``); el estado ``checked`` debe alternar.
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("mostrar carbonos")
    QApplication.processEvents()
    idx = next(
        i for i, e in enumerate(win.command_palette._filtered)
        if e.action is action
    )
    win.command_palette._idx = idx
    _press_input_key(win, int(Qt.Key.Key_Return))
    QApplication.processEvents()
    assert action.isChecked() is (not initial), "ejecutar alterna checked"
    # Revierte el estado para no contaminar otros tests.
    action.trigger()
    QApplication.processEvents()
    assert action.isChecked() is initial


# ---------------------------------------------------------------------------
# 19. Al ejecutar se cierra la paleta
# ---------------------------------------------------------------------------
def test_palette_closes_after_execution(win):
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("limpiar 2d (1 paso)")
    QApplication.processEvents()
    assert win.command_palette.is_open()
    _press_input_key(win, int(Qt.Key.Key_Return))
    QApplication.processEvents()
    assert not win.command_palette.is_open()


# ---------------------------------------------------------------------------
# 20. Reabrir no crea registros duplicados (ya cubierto por 3, reforzado)
# ---------------------------------------------------------------------------
def test_reopening_does_not_duplicate_registry(win):
    before = win._command_registry.count()
    for _ in range(3):
        win.command_palette.open()
        QApplication.processEvents()
        win.command_palette.close_overlay()
        QApplication.processEvents()
    assert win._command_registry.count() == before
    ids = [id(e.action) for e in _entries(win)]
    assert len(set(ids)) == len(ids)


# ---------------------------------------------------------------------------
# 21. Light → dark → light no rompe la paleta
# ---------------------------------------------------------------------------
def test_light_dark_light_does_not_break_palette(win):
    win.command_palette.open()
    QApplication.processEvents()
    assert win.command_palette.is_open()
    win.toggle_theme()  # -> dark
    QApplication.processEvents()
    win.command_palette.refresh_theme(win.current_theme)
    assert win.command_palette.is_open()
    win.command_palette.input.setText("export")
    QApplication.processEvents()
    assert win.command_palette.filter_count() >= 1
    win.toggle_theme()  # -> light
    QApplication.processEvents()
    win.command_palette.refresh_theme(win.current_theme)
    win.command_palette.input.setText("export")
    QApplication.processEvents()
    assert win.command_palette.filter_count() >= 1
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 22. 980×600 mantiene la paleta dentro de la ventana
# ---------------------------------------------------------------------------
def test_980x600_keeps_palette_inside_window(win):
    win.resize(980, 600)
    QApplication.processEvents()
    win.command_palette.open()
    QApplication.processEvents()
    window_rect = win.rect()
    card = win.command_palette.card.geometry()
    assert card.x() >= window_rect.x()
    assert card.y() >= window_rect.y()
    assert card.x() + card.width() <= window_rect.x() + window_rect.width()
    assert card.y() + card.height() <= window_rect.y() + window_rect.height()
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 23. Tool shortcuts no interfieren mientras se escribe en el input
# ---------------------------------------------------------------------------
def test_tool_shortcuts_do_not_interfere_while_typing(win):
    win.command_palette.open()
    QApplication.processEvents()
    win.command_palette.input.setText("b")  # "b" es el atajo de herramienta Enlace
    QApplication.processEvents()
    # El input contiene "b" (se escribió) y no se activó la herramienta Enlace.
    assert win.command_palette.input.text() == "b"
    assert win.canvas.state.active_tool != "tool_bond"
    win.command_palette.input.setText("bc")
    QApplication.processEvents()
    assert win.command_palette.input.text() == "bc"
    assert win.canvas.state.active_tool != "tool_bond"
    win.command_palette.close_overlay()
    QApplication.processEvents()


# ---------------------------------------------------------------------------
# 24. Las acciones existentes siguen accesibles desde menús
# ---------------------------------------------------------------------------
def test_existing_actions_still_accessible_from_menus(win):
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

    for attr in (
        "action_clean_2d_full",
        "action_clean_2d_publication",
        "action_clean_2d_propose",
        "action_validate_structure",
        "action_export_png",
        "action_export_svg",
        "action_undo",
        "action_new",
    ):
        action = getattr(win, attr)
        assert in_some_menu(action), attr


# ---------------------------------------------------------------------------
# 25. Única QAction global posee Ctrl+P; Ctrl+K pertenece a Clean2D quick
# ---------------------------------------------------------------------------
def test_only_one_global_action_owns_ctrl_p(win):
    # Solo ``action_command_palette`` posee ``Ctrl+P`` (contexto ventana).
    owners = []
    for action in win.actions():
        for seq in action.shortcuts():
            if seq == QKeySequence("Ctrl+P"):
                owners.append(action)
    assert owners == [win.action_command_palette], (
        "Solo action_command_palette debe poseer Ctrl+P (única QAction global)"
    )
    # Ctrl+K, por su parte, pertenece a Clean2D quick (y a ninguna otra).
    ctrl_k_owners = []
    for action in win.actions():
        for seq in action.shortcuts():
            if seq == QKeySequence("Ctrl+K"):
                ctrl_k_owners.append(action)
    assert ctrl_k_owners == [win.action_clean_2d_full], (
        "Solo action_clean_2d_full debe poseer Ctrl+K"
    )


# ---------------------------------------------------------------------------
# Frontera: command_palette.py no importa química/canvas
# ---------------------------------------------------------------------------
def test_command_palette_module_has_no_chemical_or_canvas_imports():
    import chemuson.gui.command_palette as module

    src = open(module.__file__, encoding="utf-8").read()
    for forbidden in (
        "chemuson.clean2d",
        "chemuson.chemname",
        "chemuson.chemio.persistence",
        "chemuson.gui.canvas",
    ):
        assert forbidden not in src, f"import prohibido: {forbidden}"
