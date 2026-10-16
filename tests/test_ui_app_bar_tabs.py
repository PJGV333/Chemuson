"""Tests de la barra de aplicación y las pestañas de documento (Fase 3).

OpenSpec: ``2026-09-24-modernize-ui-app-bar-document-tabs``.

Cubre:
- ``AppBar``/``SearchPill`` unitarios: presencia de controles, altura fija,
  botones que sostienen las QAction existentes;
- ``DocumentTabBar`` unitaria: espejo (sync_tabs/set_tab/select), punto de
  suciedad, botón cerrar y ``+``;
- ventana real (offscreen): pestaña inicial, documento nuevo por ``+``
  (misma ``action_new``), suciedad desde el estado real (``QUndoStack``),
  cierre/cambio/reordenamiento por el flujo existente, estado enabled de
  undo/redo, tema light→dark→light, resize 1440×900 / 980×600 y ausencia de
  atajo propio en la píldora de búsqueda (Fase 6: la píldora y Ctrl+P comparten
  el camino de apertura de la paleta de comandos; Ctrl+K volvió a Clean2D).
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PyQt6.QtGui import QUndoCommand
from PyQt6.QtWidgets import QApplication

from chemuson.gui.app_bar import AppBar, SearchPill
from chemuson.gui.document_tabs import DocumentTabBar
from chemuson.gui.theme import METRICS
from chemuson.gui.theme.icon_provider import DEFAULT_ICONS_DIR, IconProvider

REPO_ROOT = Path(__file__).resolve().parent.parent
APP_BAR_SOURCE = REPO_ROOT / "src" / "chemuson" / "gui" / "app_bar.py"
DOC_TABS_SOURCE = REPO_ROOT / "src" / "chemuson" / "gui" / "document_tabs.py"

APP_BAR_SVG_NAMES = (
    "plus",
    "search",
    "moon",
    "sun",
    "sliders",
    "x",
    "doc",
)


class _NoopCommand(QUndoCommand):
    """Comando undo mínimo para ensuciar un canvas real en tests."""

    def redo(self) -> None:
        return

    def undo(self) -> None:
        return


def _make_action_window():
    """Ventana principal real (offscreen) con la app bar montada."""
    from chemuson.gui.main_window import ChemusonWindow

    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    return win


# ---------------------------------------------------------------------------
# 1. SVG nuevos de la app bar
# ---------------------------------------------------------------------------

class TestAppBarSvgAssets:
    def test_app_bar_svgs_exist_and_are_valid(self) -> None:
        for name in APP_BAR_SVG_NAMES:
            path = DEFAULT_ICONS_DIR / f"i-{name}.svg"
            assert path.exists(), name
            root = ET.fromstring(path.read_text())
            assert root.tag.endswith("svg"), name
            assert root.get("viewBox") == "0 0 24 24", name
            assert "currentColor" in path.read_text(), name

    def test_app_bar_svgs_resolve_through_provider(self) -> None:
        provider = IconProvider()
        for name in APP_BAR_SVG_NAMES:
            for size in (14, 16, 18, 24):
                icon = provider.icon(name, "#475569", size)
                pix = icon.pixmap(size, size)
                assert not pix.isNull() and pix.width() > 0, (name, size)

    def test_app_bar_sources_only_use_known_svg_names(self) -> None:
        """Paridad: todo ``icon("name")``/``pixmap("name")`` de la app bar
        resuelve a un SVG del set (evita assets huérfanos)."""
        literals: set[str] = set()
        for source in (APP_BAR_SOURCE, DOC_TABS_SOURCE):
            literals.update(
                re.findall(r'(?s)\.(?:icon|pixmap)\(\s*"([^"]+)"', source.read_text())
            )
        for name in literals:
            assert (DEFAULT_ICONS_DIR / f"i-{name}.svg").exists(), name


# ---------------------------------------------------------------------------
# 2. Componentes unitarios (sin ventana completa)
# ---------------------------------------------------------------------------

class _FakeActions:
    """Acciones Q minimal para montar un AppBar aislado."""

    def __init__(self):
        from PyQt6.QtGui import QAction
        from PyQt6.QtWidgets import QApplication

        app = QApplication.instance()
        self.new = QAction("Nuevo", app)
        self.new.setShortcut("Ctrl+N")
        self.open = QAction("Abrir", app)
        self.open.setShortcut("Ctrl+O")
        self.save = QAction("Guardar", app)
        self.save.setShortcut("Ctrl+S")
        self.paste = QAction("Pegar", app)
        self.undo = QAction("Deshacer", app)
        self.undo.setShortcut("Ctrl+Z")
        self.redo = QAction("Rehacer", app)
        self.redo.setShortcut("Ctrl+Shift+Z")
        self.prefs = QAction("Preferencias...", app)
        self.theme = QAction("Tema", app)
        self.theme.setCheckable(True)


def _make_app_bar(actions: _FakeActions) -> AppBar:
    return AppBar(
        new_action=actions.new,
        open_action=actions.open,
        save_action=actions.save,
        paste_action=actions.paste,
        undo_action=actions.undo,
        redo_action=actions.redo,
        preferences_action=actions.prefs,
        theme_action=actions.theme,
    )


class TestAppBarUnit:
    def test_structure_and_fixed_height(self) -> None:
        app = _FakeActions()
        bar = _make_app_bar(app)
        assert bar.height() == 54 or bar.minimumHeight() == METRICS["appbarH"]
        assert bar.tab_bar is not None
        assert isinstance(bar.search_pill, SearchPill)
        assert not hasattr(bar, "brand_label")
        assert not hasattr(bar, "brand_name")
        assert not hasattr(bar, "version_label")
        assert set(bar.file_action_buttons) == {"new", "open", "save", "paste"}
        assert bar.tab_bar.new_button.isHidden()

    def test_buttons_hold_existing_qactions(self) -> None:
        app = _FakeActions()
        bar = _make_app_bar(app)
        assert bar.undo_button.defaultAction() is app.undo
        assert bar.redo_button.defaultAction() is app.redo
        assert bar.preferences_button.defaultAction() is app.prefs
        assert bar.theme_button.defaultAction() is app.theme
        assert bar.file_action_buttons["new"].defaultAction() is app.new
        assert bar.file_action_buttons["open"].defaultAction() is app.open
        assert bar.file_action_buttons["save"].defaultAction() is app.save
        assert bar.file_action_buttons["paste"].defaultAction() is app.paste

    def test_search_pill_is_entry_to_command_palette(self) -> None:
        app = _FakeActions()
        bar = _make_app_bar(app)
        # La píldora no registra QShortcut/atajo alguno: abre la command
        # palette por la misma ``QAction`` de apertura (Ctrl+P) que la ventana
        # conecta a ``searchActivated`` (un único camino de apertura).
        from PyQt6.QtGui import QShortcut

        assert bar.findChildren(QShortcut) == []
        # Entrada real (Fase 6): badge ``Ctrl P`` visible y el tooltip deja de
        # decir "(próximamente)". La píldora no es un editor permanente.
        assert bar.search_pill.text_label.text() == "Buscar o ejecutar…"
        assert "Ctrl P" in bar.search_pill.kbd.text()
        assert bar.search_pill.kbd.isVisibleTo(bar.search_pill)
        assert "(próximamente)" not in bar.search_pill.toolTip()
        from PyQt6.QtWidgets import QLineEdit

        assert bar.search_pill.findChildren(QLineEdit) == []

    def test_theme_refresh_recolors_icons(self) -> None:
        app = _FakeActions()
        bar = _make_app_bar(app)
        bar.refresh_icons("light")
        light_theme_icon = bar.theme_button.icon()
        assert not light_theme_icon.isNull()
        bar.refresh_icons("dark")
        dark_theme_icon = bar.theme_button.icon()
        assert not dark_theme_icon.isNull()
        # El botón de tema cambia de glifo (moon -> sun) sin cambiar la altura.
        assert light_theme_icon.cacheKey() != dark_theme_icon.cacheKey()
        bar.refresh_icons("light")
        assert bar.theme_button.icon().cacheKey() == light_theme_icon.cacheKey()


class TestDocumentTabBarUnit:
    def _bar(self) -> DocumentTabBar:
        return DocumentTabBar("#475569")

    def test_sync_tabs_builds_mirror(self) -> None:
        bar = self._bar()
        bar.sync_tabs(["A.cmsn", "Sin título 2"], [True, False], 1)
        assert bar.tab_count() == 2
        assert bar.tab_bar.tabText(0) == "A.cmsn"
        assert bar.tab_bar.tabText(1) == "Sin título 2"
        assert bar.current_index() == 1
        assert bar._side_button(0).is_dirty_shown()
        assert not bar._side_button(1).is_dirty_shown()

    def test_sync_tabs_is_idempotent(self) -> None:
        bar = self._bar()
        bar.sync_tabs(["A"], [False], 0)
        bar.sync_tabs(["A", "B"], [True, False], 0)
        assert bar.tab_count() == 2
        assert bar.tab_bar.tabText(0) == "A"

    def test_set_tab_updates_dirty_dot(self) -> None:
        bar = self._bar()
        bar.sync_tabs(["A"], [False], 0)
        bar.set_tab(0, "A", True)
        assert bar._side_button(0).is_dirty_shown()
        bar.set_tab(0, "A", False)
        assert not bar._side_button(0).is_dirty_shown()

    def test_close_button_emits_current_index(self) -> None:
        bar = self._bar()
        bar.sync_tabs(["A", "B", "C"], [False, False, False], 0)
        bar.show()
        QApplication.processEvents()
        got: list[int] = []
        bar.closeRequested.connect(got.append)
        side = bar._side_button(1)
        side.close_button.click()
        assert got == [1]

    def test_new_button_emits_signal(self) -> None:
        bar = self._bar()
        fired: list[bool] = []
        bar.newDocumentRequested.connect(lambda: fired.append(True))
        bar.new_button.click()
        assert fired == [True]

    def test_tab_click_emits_tab_activated(self) -> None:
        bar = self._bar()
        bar.sync_tabs(["A", "B"], [False, False], 0)
        got: list[int] = []
        bar.tabActivated.connect(got.append)
        bar.tab_bar.tabBarClicked.emit(1)
        assert got == [1]


# ---------------------------------------------------------------------------
# 3. Ventana real
# ---------------------------------------------------------------------------

class TestRealWindow:
    def test_app_bar_mounted_and_main_toolbar_hidden_by_default(
        self, isolated_shell_settings
    ) -> None:
        win = _make_action_window()
        assert win.app_bar is not None and win.app_bar.isVisibleTo(win)
        assert win.app_bar.minimumHeight() == 54
        assert win.tabs.tabBar().isHidden()
        assert win.main_toolbar.isHidden()
        assert win.action_show_menu_bar.isChecked()
        assert win.action_show_menu_bar in win.view_menu.actions()
        # La AppBar integra cuatro acciones sin repetir controladores.
        assert win.app_bar.file_action_buttons["new"].defaultAction() is win.action_new
        assert win.app_bar.file_action_buttons["open"].defaultAction() is win.action_open
        assert win.app_bar.file_action_buttons["save"].defaultAction() is win.action_save
        assert win.app_bar.file_action_buttons["paste"].defaultAction() is win.action_paste
        assert win.app_bar.tab_bar.new_button.isHidden()
        assert win.main_toolbar.isHidden()
        # La barra opcional conserva sus tres acciones históricas.
        assert [a for a in win.main_toolbar.actions() if not a.isSeparator()] == [
            win.action_new,
            win.action_open,
            win.action_save,
        ]
        assert win.action_show_menu_bar.text() == "Mostrar barra de herramientas"
        assert win.view_menu.actions().count(win.action_show_menu_bar) == 1
        assert win.main_toolbar.toggleViewAction().isEnabled() is False
        # La QMenuBar tradicional con los seis menús es visible por defecto
        # (compacta) y sigue siendo la misma: los seis QMenu también están en
        # el menú agregado de la hamburguesa.
        assert win.menuBar() is not None
        assert win.menuBar().isVisible() is True
        assert win.menuBar().height() <= 36
        menu_titles = [a.text() for a in win.menuBar().actions()]
        for title in ("Archivo", "Editar", "Ver", "Estructura", "Reacción", "Ayuda"):
            assert any(title in m for m in menu_titles), title
        # Wrapper central [app_bar, text_toolbar, body(rail, splitter(canvas, side_panel))]
        central = win.centralWidget()
        assert central is not win.tabs
        assert win.app_bar.parent() is central
        assert win._central_splitter.parent() is central
        assert win.tabs.parent() is win._central_splitter
        assert win.side_panel.parent() is win._central_splitter
        assert win.side_panel.isHidden()
        assert win.text_toolbar.parent() is central
        assert win.tool_rail.parent() is central

    def test_vertical_hierarchy_appbar_above_text_toolbar(
        self, isolated_shell_settings
    ) -> None:
        """Convergencia con el spike: AppBar -> TextFormatToolbar (oculta
        por defecto) -> [rail | splitter(canvas | side panel)].

        La ``TextFormatToolbar`` ya no es toolbar de área superior de
        ``QMainWindow``; vive en el layout central, debajo de la app bar y
        sobre la franja canvas, conservando sus acciones y señales. La
        ``QMenuBar`` tradicional es visible por defecto y, si se oculta,
        los menús siguen accesibles por la hamburguesa/Alt.
        """
        win = _make_action_window()
        assert win.menuBar().isVisible() is True
        # La toolbar de texto es contextual: oculta por defecto.
        assert win.text_toolbar.isVisible() is False
        central = win.centralWidget()
        layout = central.layout()
        assert layout is not None
        widgets = [
            layout.itemAt(i).widget()
            for i in range(layout.count())
            if layout.itemAt(i) is not None and layout.itemAt(i).widget() is not None
        ]
        # Orden vertical exacto dentro del wrapper central:
        # app_bar → text_toolbar (oculta) → banda body (rail + splitter).
        assert win.app_bar in widgets
        assert win.text_toolbar in widgets
        assert widgets.index(win.app_bar) < widgets.index(win.text_toolbar)
        body_item = layout.itemAt(2)
        assert body_item is not None and body_item.layout() is not None
        body = body_item.layout()
        body_widgets = [
            body.itemAt(i).widget()
            for i in range(body.count())
            if body.itemAt(i) is not None and body.itemAt(i).widget() is not None
        ]
        assert win.tool_rail in body_widgets
        assert win._central_splitter in body_widgets
        assert body_widgets.index(win.tool_rail) < body_widgets.index(win._central_splitter)
        assert win._central_splitter.widget(0) is win.tabs
        assert win._central_splitter.widget(1) is win.side_panel
        # No es una toolbar de área de QMainWindow (no compite por la franja
        # superior con la app bar): su padre es el wrapper central, no la
        # ventana.
        from PyQt6.QtWidgets import QToolBar

        direct_toolbars = [
            tb for tb in win.findChildren(QToolBar) if tb.parentWidget() is win
        ]
        assert win.text_toolbar not in direct_toolbars
        assert win.toolbar in direct_toolbars  # la izquierda sigue montada
        # El handler real sigue en la ventana y la señal sigue viva.
        assert callable(win._on_text_format_changed)
        fired: list[bool] = []
        win.text_toolbar.format_changed.connect(lambda *args: fired.append(True))
        win.text_toolbar.format_changed.emit("Arial", 12, True, False, False, False, False, "bold")
        assert fired == [True]

    def test_initial_single_tab_clean(self) -> None:
        win = _make_action_window()
        assert win.tabs.count() == 1
        assert win.app_bar.tab_count() == 1
        assert win.app_bar.tab_bar.tab_bar.tabText(0) == "Sin título"
        assert not win.app_bar.tab_bar._side_button(0).is_dirty_shown()
        assert win.app_bar.tab_bar.current_index() == 0

    def test_new_document_via_appbar_action_button(self) -> None:
        win = _make_action_window()
        fired = []
        # El botón Nuevo de la AppBar reutiliza la misma QAction.
        win.action_new.triggered.connect(lambda: fired.append(True))
        win.app_bar.file_action_buttons["new"].click()
        QApplication.processEvents()
        assert fired == [True]
        assert win.tabs.count() == 2
        assert win.app_bar.tab_count() == 2
        assert win.app_bar.tab_bar.tab_bar.tabText(1) == "Sin título 2"
        assert win.app_bar.tab_bar.current_index() == 1
        assert win.canvas is win._canvas_from_tab_index(1)

    def test_dirty_state_from_real_undo_stack(self) -> None:
        win = _make_action_window()
        bar = win.app_bar.tab_bar
        assert not bar._side_button(0).is_dirty_shown()
        win.canvas.undo_stack.push(_NoopCommand())
        QApplication.processEvents()
        # Contrato histórico: el QTabWidget añade " *" (update_tab_title)...
        assert win.tabs.tabText(0) == "Sin título *"
        # ...y el espejo muestra el punto de suciedad (mismo estado real).
        assert bar._side_button(0).is_dirty_shown()
        # Limpieza real -> ambos vuelven.
        win.canvas.undo_stack.setClean()
        QApplication.processEvents()
        assert win.tabs.tabText(0) == "Sin título"
        assert not bar._side_button(0).is_dirty_shown()

    def test_close_and_switch_tabs(self) -> None:
        win = _make_action_window()
        win._on_file_new()
        QApplication.processEvents()
        assert win.tabs.count() == 2
        # Cambiar por la pestaña 1 del espejo (flujo existente).
        win.app_bar.tab_bar.tab_bar.tabBarClicked.emit(0)
        QApplication.processEvents()
        assert win.canvas is win._canvas_from_tab_index(0)
        # Cerrar la pestaña 1 (limpia: sin diálogo).
        win.app_bar.closeRequested.emit(1)
        QApplication.processEvents()
        assert win.tabs.count() == 1
        assert win.app_bar.tab_count() == 1
        assert win.canvas is win._canvas_from_tab_index(0)

    def test_reorder_tabs(self) -> None:
        win = _make_action_window()
        win._on_file_new()
        win._on_file_new()
        QApplication.processEvents()
        assert win.tabs.count() == 3
        titles_before = [win.tabs.tabText(i) for i in range(3)]
        assert titles_before[0] != titles_before[2]
        # Drag en el espejo -> moveTab en el QTabWidget (mismo contrato).
        win.app_bar.tabMoved.emit(0, 2)
        QApplication.processEvents()
        titles_after = [win.tabs.tabText(i) for i in range(3)]
        assert titles_after[2] == titles_before[0]
        assert win.app_bar.tab_count() == 3

    def test_shared_qactions_and_enabled_states(self) -> None:
        win = _make_action_window()
        bar = win.app_bar
        assert bar.undo_button.defaultAction() is win.action_undo
        assert bar.redo_button.defaultAction() is win.action_redo
        assert bar.preferences_button.defaultAction() is win.action_preferences
        assert bar.file_action_buttons["new"].defaultAction() is win.action_new
        assert bar.file_action_buttons["open"].defaultAction() is win.action_open
        assert bar.file_action_buttons["save"].defaultAction() is win.action_save
        assert bar.file_action_buttons["paste"].defaultAction() is win.action_paste
        assert win.action_paste not in win.main_toolbar.actions()
        edit_menu = next(
            action.menu()
            for action in win.menuBar().actions()
            if action.text() == "Editar"
        )
        assert edit_menu.actions().count(win.action_paste) == 1
        # Atajos originales intactos, sin duplicados (estándar de plataforma).
        from PyQt6.QtGui import QKeySequence

        assert win.action_undo.shortcut().toString() == QKeySequence(
            QKeySequence.StandardKey.Undo
        ).toString()
        assert win.action_redo.shortcut().toString() == QKeySequence(
            QKeySequence.StandardKey.Redo
        ).toString()
        # Estado real: sin historial, undo/redo deshabilitados en botón y acción.
        assert not win.action_undo.isEnabled()
        assert not bar.undo_button.isEnabled()
        assert not win.action_redo.isEnabled()
        assert not bar.redo_button.isEnabled()
        # Un comando real habilita la acción y, con ella, el botón.
        win.canvas.undo_stack.push(_NoopCommand())
        QApplication.processEvents()
        assert win.action_undo.isEnabled()
        assert bar.undo_button.isEnabled()
        win.canvas.undo_stack.setClean()

    def test_search_pill_has_no_shortcut_and_shortcuts_split(self) -> None:
        win = _make_action_window()
        from PyQt6.QtGui import QShortcut

        # La píldora no registra QShortcut propio: comparte el camino de
        # apertura (``action_command_palette``) con Ctrl+P (Fase 6).
        assert win.app_bar.findChildren(QShortcut) == []
        # Ctrl+P es el atajo global de la paleta; Ctrl+K volvió a pertenecer a
        # ``action_clean_2d_full`` (corrección de UX).
        assert win.action_command_palette.shortcut().toString() == "Ctrl+P"
        assert win.action_clean_2d_full.shortcut().toString() == "Ctrl+K"
        assert win.action_save.shortcut().toString() == "Ctrl+S"

    def test_theme_cycle_light_dark_light(self) -> None:
        win = _make_action_window()
        # Normaliza el estado inicial (el tema persistido por M21 puede
        # llegar en cualquiera de los dos valores).
        win.current_theme = "light"
        win._apply_theme()
        QApplication.processEvents()
        bar = win.app_bar
        light_theme_icon = bar.theme_button.icon().cacheKey()
        win.toggle_theme()
        QApplication.processEvents()
        assert win.action_theme_toggle.isChecked()
        assert bar.theme_button.icon().cacheKey() != light_theme_icon
        win.toggle_theme()
        QApplication.processEvents()
        assert not win.action_theme_toggle.isChecked()
        # Vuelta exacta a claro (caché theme-aware, sin contaminación).
        assert bar.theme_button.icon().cacheKey() == light_theme_icon
        # La altura no cambia con el tema.
        assert bar.minimumHeight() == 54

    def test_resize_targets(self) -> None:
        win = _make_action_window()
        for _ in range(3):
            win._on_file_new()
        QApplication.processEvents()
        win.resize(1440, 900)
        QApplication.processEvents()
        assert win.app_bar.minimumHeight() == 54
        assert win.app_bar.tab_count() == 4
        win.resize(980, 600)
        QApplication.processEvents()
        assert win.app_bar.minimumHeight() == 54
        assert win.app_bar.tab_count() == 4
        assert win.tabs.isVisible()
        # El canvas activo sigue montado.
        assert win.canvas is win._canvas_from_tab_index(win.tabs.currentIndex())
