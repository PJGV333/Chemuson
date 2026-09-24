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
  atajo en la píldora de búsqueda (Ctrl+K sigue siendo "Clean 2D full").
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
    "flask",
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


def _pixmap_signature(pixmap) -> list[int]:
    """Firma de píxeles para comparar QPixmaps (PyQt6 no define ``==``)."""
    image = pixmap.toImage()
    w, h = image.width(), image.height()
    return [
        image.pixelColor(x, y).rgba()
        for x, y in ((0, 0), (w // 2, h // 2), (max(w - 1, 0), max(h - 1, 0)), (3, 8), (w - 4, 5))
    ]


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
        self.undo = QAction("Deshacer", app)
        self.undo.setShortcut("Ctrl+Z")
        self.redo = QAction("Rehacer", app)
        self.redo.setShortcut("Ctrl+Shift+Z")
        self.prefs = QAction("Preferencias...", app)
        self.theme = QAction("Tema", app)
        self.theme.setCheckable(True)


class TestAppBarUnit:
    def test_structure_and_fixed_height(self) -> None:
        app = _FakeActions()
        bar = AppBar(
            version="0.4",
            undo_action=app.undo,
            redo_action=app.redo,
            preferences_action=app.prefs,
            theme_action=app.theme,
        )
        assert bar.height() == 54 or bar.minimumHeight() == METRICS["appbarH"]
        assert bar.tab_bar is not None
        assert isinstance(bar.search_pill, SearchPill)
        assert not bar.brand_label.pixmap().isNull()
        assert bar.version_label.text() == "0.4"

    def test_buttons_hold_existing_qactions(self) -> None:
        app = _FakeActions()
        bar = AppBar(
            version="0.4",
            undo_action=app.undo,
            redo_action=app.redo,
            preferences_action=app.prefs,
            theme_action=app.theme,
        )
        assert bar.undo_button.defaultAction() is app.undo
        assert bar.redo_button.defaultAction() is app.redo
        assert bar.preferences_button.defaultAction() is app.prefs
        assert bar.theme_button.defaultAction() is app.theme

    def test_search_pill_is_placeholder_without_shortcut(self) -> None:
        app = _FakeActions()
        bar = AppBar(
            version="0.4",
            undo_action=app.undo,
            redo_action=app.redo,
            preferences_action=app.prefs,
            theme_action=app.theme,
        )
        # La píldora no registra QShortcut/atajo alguno (Ctrl+K es de
        # ``action_clean_2d_full``; la command palette es la Fase 6).
        from PyQt6.QtGui import QShortcut

        assert bar.findChildren(QShortcut) == []
        assert "Ctrl K" in bar.search_pill.kbd.text()

    def test_theme_refresh_recolors_icons(self) -> None:
        app = _FakeActions()
        bar = AppBar(
            version="0.4",
            undo_action=app.undo,
            redo_action=app.redo,
            preferences_action=app.prefs,
            theme_action=app.theme,
        )
        bar.refresh_icons("light")
        light_theme_icon = bar.theme_button.icon()
        light_brand = bar.brand_label.pixmap().copy()
        assert not light_theme_icon.isNull()
        bar.refresh_icons("dark")
        dark_theme_icon = bar.theme_button.icon()
        assert not dark_theme_icon.isNull()
        # El botón de tema cambia de glifo (moon -> sun) y el tinte de la
        # marca pasa a usar el accent oscuro.
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
    def test_app_bar_mounted_and_classic_toolbar_hidden(self) -> None:
        win = _make_action_window()
        assert win.app_bar is not None and win.app_bar.isVisibleTo(win)
        assert win.app_bar.minimumHeight() == 54
        assert win.tabs.tabBar().isHidden()
        assert not win.main_toolbar.isVisible()
        assert win.menuBar().isVisible()
        # Los 6 menús históricos siguen presentes.
        menu_titles = [a.text() for a in win.menuBar().actions()]
        for title in ("Archivo", "Editar", "Ver", "Estructura", "Reacción", "Ayuda"):
            assert any(title in m for m in menu_titles), title
        # Wrapper central [app_bar, tabs]
        central = win.centralWidget()
        assert central is not win.tabs
        assert win.app_bar.parent() is central
        assert win.tabs.parent() is central

    def test_initial_single_tab_clean(self) -> None:
        win = _make_action_window()
        assert win.tabs.count() == 1
        assert win.app_bar.tab_count() == 1
        assert win.app_bar.tab_bar.tab_bar.tabText(0) == "Sin título"
        assert not win.app_bar.tab_bar._side_button(0).is_dirty_shown()
        assert win.app_bar.tab_bar.current_index() == 0

    def test_new_document_via_plus_button(self) -> None:
        win = _make_action_window()
        fired = []
        # La app bar reutiliza la misma action_new (no un segundo handler).
        win.action_new.triggered.connect(lambda: fired.append(True))
        win.app_bar.tab_bar.new_button.click()
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

    def test_search_pill_has_no_shortcut_and_ctrl_k_unchanged(self) -> None:
        win = _make_action_window()
        from PyQt6.QtGui import QShortcut

        assert win.app_bar.findChildren(QShortcut) == []
        # Ctrl+K sigue siendo "Clean 2D full" (la píldora no lo reclama).
        assert win.action_clean_2d_full.shortcut().toString() == "Ctrl+K"

    def test_theme_cycle_light_dark_light(self) -> None:
        win = _make_action_window()
        # Normaliza el estado inicial (el tema persistido por M21 puede
        # llegar en cualquiera de los dos valores).
        win.current_theme = "light"
        win._apply_theme()
        QApplication.processEvents()
        bar = win.app_bar
        light_theme_icon = bar.theme_button.icon().cacheKey()
        light_brand = bar.brand_label.pixmap().copy()
        win.toggle_theme()
        QApplication.processEvents()
        assert win.action_theme_toggle.isChecked()
        assert bar.theme_button.icon().cacheKey() != light_theme_icon
        assert not bar.brand_label.pixmap().isNull()
        win.toggle_theme()
        QApplication.processEvents()
        assert not win.action_theme_toggle.isChecked()
        # Vuelta exacta a claro (caché theme-aware, sin contaminación).
        assert bar.theme_button.icon().cacheKey() == light_theme_icon
        assert _pixmap_signature(bar.brand_label.pixmap()) == _pixmap_signature(light_brand)
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
