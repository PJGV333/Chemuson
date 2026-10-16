"""Correcciones de shell de la beta: AppBar de dos franjas y pestañas estables.

OpenSpec: ``fix-beta-tab-close-and-toolbar-visibility`` (Fase A de la
corrección solicitada tras la aceptación visual).

Cubre:
- A1: la barra de menús tradicional (Archivo/Editar/Ver/Estructura/Reacción/
  Ayuda) es visible por defecto y compacta, reutiliza los mismos ``QMenu`` y
  ``QAction``, puede ocultarse desde Ver conservando el acceso hamburguesa, y
  el default nuevo migra una sola vez sobre preferencias previas;
- A2: «Mostrar barra de herramientas» controla sólo la QMenuBar; la toolbar
  histórica permanece oculta incluso con preferencias antiguas;
- A3: pasar de pestaña limpia a modificada no altera la geometría de la
  pestaña ni la posición de la X; el punto dirty sólo se pinta cuando toca y
  el cierre sigue apuntando al documento correcto.
"""

from __future__ import annotations

from collections import Counter

import pytest
from PyQt6.QtCore import QPoint, QRect, QSettings, Qt
from PyQt6.QtGui import QColor, QKeySequence, QPainter, QPen, QPixmap
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QMessageBox

from chemuson.gui.document_tabs import (
    _CLOSE_BTN_SIZE,
    _DOT_SIZE,
    _SIDE_MARGIN_X,
    _SIDE_SPACING,
    _SIDE_WIDTH,
    DocumentTabBar,
)
from chemuson.gui.main_window import ChemusonWindow
from chemuson.platform.settings import (
    SHELL_DEFAULTS_VERSION,
    load_ui_preferences,
)

MENU_TITLES = ("Archivo", "Editar", "Ver", "Estructura", "Reacción", "Ayuda")
#: Techo de altura de la fila de menús (baseline medido: 43 px → objetivo ≤36).
MAX_MENU_BAR_HEIGHT = 36


def _pin_settings(tmp_path, monkeypatch, values: dict[str, object] | None = None):
    """Fija un INI temporal determinista para contratos de *defaults*."""
    settings_path = tmp_path / "Chemuson-shell.ini"
    seed = QSettings(str(settings_path), QSettings.Format.IniFormat)
    seed.clear()
    seed.setValue("ui/onboarding/completed", True)
    for key, value in (values or {}).items():
        seed.setValue(key, value)
    seed.sync()

    def factory():
        # Los QSettings del mismo archivo comparten estado: limpiar en cada
        # llamada borraría lo que la ventana acaba de persistir.
        return QSettings(str(settings_path), QSettings.Format.IniFormat)

    monkeypatch.setattr("chemuson.platform.settings.application_settings", factory)
    monkeypatch.setattr("chemuson.gui.shell.assembly.application_settings", factory)
    return factory


def _open_window(tmp_path, monkeypatch, values: dict[str, object] | None = None):
    _pin_settings(tmp_path, monkeypatch, values)
    win = ChemusonWindow()
    win.resize(1280, 900)
    win.show()
    QApplication.processEvents()
    return win


@pytest.fixture
def shell(tmp_path, monkeypatch):
    """Ventana real con menú visible y toolbar opcional oculta por defecto."""
    win = _open_window(tmp_path, monkeypatch)
    yield win
    for index in range(win.tabs.count()):
        canvas = win._canvas_from_tab_index(index)
        if canvas is not None:
            canvas.undo_stack.setClean()
    win.close()
    QApplication.processEvents()


def _close_window(win: ChemusonWindow) -> None:
    for index in range(win.tabs.count()):
        canvas = win._canvas_from_tab_index(index)
        if canvas is not None:
            canvas.undo_stack.setClean()
    win.close()
    QApplication.processEvents()


def _top_menus(win: ChemusonWindow) -> list:
    return [a.menu() for a in win.menuBar().actions() if a.menu() is not None]


def _global_widget_rect(widget, local_rect: QRect) -> QRect:
    return QRect(widget.mapToGlobal(local_rect.topLeft()), local_rect.size())


def _assert_popup_fits_available_screen(popup, available: QRect) -> None:
    popup_rect = _global_widget_rect(popup, popup.rect())
    assert available.contains(popup_rect), (
        f"popup {popup.objectName() or popup.title()!r} rect={popup_rect}, "
        f"available={available}"
    )
    assert popup.height() <= available.height()
    actions = [
        action
        for action in popup.actions()
        if action.isVisible() and not action.isSeparator()
    ]
    assert actions
    for action in (actions[0], actions[-1]):
        row = _global_widget_rect(popup, popup.actionGeometry(action))
        assert available.contains(row), (
            f"menu row {action.text()!r} rect={row}, available={available}"
        )


def _save_popup_screen_evidence(popup, screen, destination) -> None:
    available = screen.availableGeometry()
    image = QPixmap(available.size())
    image.fill(QColor("#f2f3f5"))
    painter = QPainter(image)
    painter.setPen(QPen(QColor("#667085"), 2))
    painter.drawRect(image.rect().adjusted(1, 1, -2, -2))
    origin = popup.mapToGlobal(popup.rect().topLeft()) - available.topLeft()
    painter.drawPixmap(origin, popup.grab())
    painter.end()
    assert image.save(str(destination))


def _mark_dirty(canvas, element: str = "C") -> None:
    canvas.model.add_atom(element, 12.0, 18.0)
    canvas._rebuild_items_from_model()
    canvas.undo_stack.resetClean()
    QApplication.processEvents()


def _rendered_colors(widget) -> Counter[str]:
    image = widget.grab().toImage()
    return Counter(
        image.pixelColor(x, y).name()
        for x in range(image.width())
        for y in range(image.height())
        if image.pixelColor(x, y).alpha() > 0
    )


def _side_snapshot(bar, index: int) -> tuple:
    """Geometría observable de la pestaña: rect, bloque lateral y X.

    El último elemento es la holgura entre el borde derecho de la X y el
    borde derecho de la pestaña: es la magnitud que debe permanecer fija
    (la posición absoluta cambia si cambia el ancho del texto, y eso no es
    un desplazamiento de la X dentro de su pestaña).
    """
    side = bar._side_button(index)
    assert side is not None
    tab_rect = bar.tab_bar.tabRect(index)
    close_in_tab = side.mapTo(bar.tab_bar, side.close_button.pos())
    right_gap = (
        tab_rect.right() - (close_in_tab.x() + side.close_button.width()) + 1
    )
    return (
        tab_rect,
        side.size(),
        side.close_button.size(),
        side.close_button.pos(),
        close_in_tab,
        right_gap,
    )


def _close_fits_in_tab(snapshot: tuple) -> bool:
    """La X completa debe quedar dentro del rect de su pestaña."""
    tab_rect, _, close_size, _, close_in_tab, _ = snapshot
    return tab_rect.contains(
        QRect(
            close_in_tab.x(),
            close_in_tab.y(),
            close_size.width(),
            close_size.height(),
        )
    )


def _invariants(snapshot: tuple) -> tuple:
    """Parte de la geometría que no debe cambiar nunca.

    Excluye la posición absoluta de la X dentro de la barra: esa sí cambia
    cuando cambia el ancho del título (elisión) o el ancho disponible, y no
    es un desplazamiento de la X dentro de su pestaña.
    """
    _, side_size, close_size, close_pos, _, right_gap = snapshot
    return (side_size, close_size, close_pos, right_gap)


def _documents(win: ChemusonWindow) -> list:
    return [win._canvas_from_tab_index(i) for i in range(win.tabs.count())]


# ---------------------------------------------------------------------------
# A1 — barra de menús tradicional
# ---------------------------------------------------------------------------


def test_traditional_menu_bar_is_visible_by_default_and_compact(shell) -> None:
    menubar = shell.menuBar()
    assert shell.action_show_menu_bar.isChecked()
    assert menubar.isVisible()
    assert [menu.title() for menu in _top_menus(shell)] == list(MENU_TITLES)
    # Fila compacta: la barra de menús no duplica la altura de la app bar.
    assert menubar.height() <= MAX_MENU_BAR_HEIGHT
    assert shell.app_bar.height() == 54
    # Sin solapamientos: la barra de menús termina antes de la app bar.
    menubar_bottom = menubar.mapTo(shell, QPoint(0, menubar.height())).y()
    app_bar_top = shell.app_bar.mapTo(shell, QPoint(0, 0)).y()
    assert menubar_bottom <= app_bar_top
    # Cada menú ocupa un hueco real y legible dentro de la barra.
    for action in menubar.actions():
        geometry = menubar.actionGeometry(action)
        assert geometry.width() > 20 and geometry.height() > 12
        assert menubar.rect().contains(geometry.topLeft())


def test_menu_bar_reuses_the_same_menus_and_actions_without_duplicates(shell) -> None:
    menus = _top_menus(shell)
    aggregated = [
        action.menu()
        for action in shell._menu_aggregator.actions()
        if action.menu() is not None
    ]
    # La hamburguesa agrega los mismos QMenu: no hay una segunda copia.
    assert aggregated == menus
    assert shell.action_save in menus[0].actions()
    assert shell.action_undo in menus[1].actions()
    assert shell.action_clean_2d in menus[3].actions()
    # Una sola instancia de cada QAction en la barra de menús.
    flat = [action for menu in menus for action in menu.actions()]
    duplicated = [action for action in set(flat) if flat.count(action) > 1]
    assert not duplicated


def test_view_popup_and_submenus_fit_available_screen(shell, tmp_path) -> None:
    menu = shell.view_menu
    screen = shell.screen()
    available = screen.availableGeometry()
    view_action = next(
        action
        for action in shell.menuBar().actions()
        if action.menu() is menu
    )

    # Reuse one real window/QMenu across the matrix. Qt's offscreen platform
    # has a fixed QScreen, so resize the app viewport and compare popup
    # coordinates to the actual availableGeometry reported by that backend.
    for theme in ("light", "dark"):
        shell.current_theme = theme
        shell._apply_theme()
        for window_size in ((1600, 900), (1280, 720), (1024, 768)):
            shell.resize(*window_size)
            QApplication.processEvents()
            QTest.mouseClick(
                shell.menuBar(),
                Qt.MouseButton.LeftButton,
                pos=shell.menuBar().actionGeometry(view_action).center(),
            )
            QApplication.processEvents()

            assert menu.isVisible()
            assert menu.height() <= available.height()
            visible = [
                action
                for action in menu.actions()
                if action.isVisible() and not action.isSeparator()
            ]
            assert visible[0] is shell.action_show_menu_bar
            assert shell.action_show_menu_bar.text() == "Mostrar barra de herramientas"
            assert shell.action_show_menu_bar.isCheckable()
            assert shell.action_show_menu_bar.isChecked()
            _assert_popup_fits_available_screen(menu, available)

            if theme == "light" and window_size == (1600, 900):
                evidence = tmp_path / "view-menu-screen-bounds.png"
                _save_popup_screen_evidence(menu, screen, evidence)
                assert evidence.stat().st_size > 0

            for parent_action in visible:
                submenu = parent_action.menu()
                if submenu is None:
                    continue
                QTest.mouseClick(
                    menu,
                    Qt.MouseButton.LeftButton,
                    pos=menu.actionGeometry(parent_action).center(),
                )
                QTest.qWait(350)
                QApplication.processEvents()
                assert submenu.isVisible(), (
                    f"theme={theme}, window_size={window_size}, "
                    f"submenu={parent_action.text()}"
                )
                assert submenu.height() <= available.height()
                _assert_popup_fits_available_screen(submenu, available)
                submenu.close()
                QApplication.processEvents()

            menu.close()
            QApplication.processEvents()


def test_view_menu_submenus_reuse_actions_and_remove_obsolete_auxiliary_option(
    shell,
) -> None:
    menus = {
        action.text(): action.menu()
        for action in shell.view_menu.actions()
        if action.menu() is not None
    }
    visible = [
        action
        for action in shell.view_menu.actions()
        if action.isVisible() and not action.isSeparator()
    ]
    assert visible[0] is shell.action_show_menu_bar
    assert [action.text() for action in visible] == [
        "Mostrar barra de herramientas",
        "Visualización",
        "Zoom",
        "Numeración",
        "Tamaño de lienzo",
        "Paletas",
        "Paneles",
    ]
    assert not hasattr(shell, "action_show_main_toolbar_aux")
    assert menus["Visualización"].actions() == [
        shell.action_show_carbons,
        shell.action_show_hydrogens,
        shell.action_aromatic_circles,
        shell.action_style,
        shell.action_rules,
        shell.action_crosshair,
    ]
    assert menus["Zoom"].actions() == [
        shell.action_zoom_in,
        shell.action_zoom_out,
        shell.action_zoom_reset,
    ]
    assert shell.action_copy in _top_menus(shell)[1].actions()
    assert shell.action_copy.shortcut() == QKeySequence(
        QKeySequence.StandardKey.Copy
    )
    assert shell.action_zoom_in in menus["Zoom"].actions()
    assert shell.action_zoom_out in menus["Zoom"].actions()
    assert shell.action_zoom_in.shortcut() == QKeySequence(
        QKeySequence.StandardKey.ZoomIn
    )
    assert shell.action_zoom_out.shortcut() == QKeySequence(
        QKeySequence.StandardKey.ZoomOut
    )
    assert shell.action_toggle_side_panel in menus["Paneles"].actions()
    assert all(
        action in menus["Paneles"].actions()
        for action in shell.side_panel_actions.values()
    )


def test_mouse_toggle_restores_menu_bar_from_hamburger_and_persists(
    tmp_path, monkeypatch
) -> None:
    settings_factory = _pin_settings(tmp_path, monkeypatch)
    win = ChemusonWindow()
    win.resize(1280, 720)
    win.show()
    QApplication.processEvents()
    menu = win.view_menu
    action = win.action_show_menu_bar
    try:
        top_action = next(
            item for item in win.menuBar().actions() if item.menu() is menu
        )
        QTest.mouseClick(
            win.menuBar(),
            Qt.MouseButton.LeftButton,
            pos=win.menuBar().actionGeometry(top_action).center(),
        )
        QApplication.processEvents()
        QTest.mouseClick(
            menu,
            Qt.MouseButton.LeftButton,
            pos=menu.actionGeometry(action).center(),
        )
        QApplication.processEvents()
        assert not win.menuBar().isVisible()
        assert not action.isChecked()
        assert load_ui_preferences(settings_factory()).show_menu_bar is False
        assert win.app_bar.isVisible()
        assert win.main_toolbar.isHidden()
        menu.close()
        win._menu_aggregator.close()
    finally:
        menu.close()
        win._menu_aggregator.close()
        _close_window(win)

    reopened = ChemusonWindow()
    reopened.resize(1280, 720)
    reopened.show()
    QApplication.processEvents()
    try:
        assert not reopened.menuBar().isVisible()
        assert not reopened.action_show_menu_bar.isChecked()
        assert reopened.main_toolbar.isHidden()

        QTest.mouseClick(reopened.app_bar.menu_button, Qt.MouseButton.LeftButton)
        QApplication.processEvents()
        aggregator = reopened._menu_aggregator
        assert aggregator.isVisible()
        view_action = next(
            item for item in aggregator.actions() if item.menu() is reopened.view_menu
        )
        QTest.mouseMove(
            aggregator,
            aggregator.actionGeometry(view_action).center(),
            delay=250,
        )
        QTest.qWait(350)
        QApplication.processEvents()
        assert reopened.view_menu.isVisible()
        restore_action = reopened.action_show_menu_bar
        QTest.mouseClick(
            reopened.view_menu,
            Qt.MouseButton.LeftButton,
            pos=reopened.view_menu.actionGeometry(restore_action).center(),
        )
        QApplication.processEvents()
        assert reopened.menuBar().isVisible()
        assert restore_action.isChecked()
        assert reopened.main_toolbar.isHidden()
        assert load_ui_preferences(reopened._settings).show_menu_bar is True
        reopened.view_menu.close()
        aggregator.close()
    finally:
        reopened.view_menu.close()
        reopened._menu_aggregator.close()
        _close_window(reopened)

    persisted = ChemusonWindow()
    persisted.show()
    QApplication.processEvents()
    try:
        assert persisted.menuBar().isVisible()
        assert persisted.action_show_menu_bar.isChecked()
        assert persisted.main_toolbar.isHidden()
    finally:
        _close_window(persisted)


def test_default_header_has_one_compact_action_row_and_shared_file_actions(shell, tmp_path) -> None:
    from PyQt6.QtWidgets import QLabel

    for width in (1440, 980, 900):
        shell.resize(width, 600)
        QApplication.processEvents()
        bar = shell.app_bar
        assert shell.menuBar().isVisible()
        assert bar.height() == 54
        assert shell.main_toolbar.isHidden()
        assert shell.action_show_menu_bar.isChecked()
        assert shell.action_show_menu_bar.text() == "Mostrar barra de herramientas"
        assert not shell.main_toolbar.toggleViewAction().isEnabled()
        assert not bar.findChild(QLabel, "appBrandName")
        assert not bar.findChild(QLabel, "appVersionPill")
        assert bar.tab_bar.new_button.isHidden()
        assert bar.file_action_buttons["new"].defaultAction() is shell.action_new
        assert bar.file_action_buttons["open"].defaultAction() is shell.action_open
        assert bar.file_action_buttons["save"].defaultAction() is shell.action_save
        assert bar.file_action_buttons["paste"].defaultAction() is shell.action_paste
        for button in bar.file_action_buttons.values():
            assert button.isVisible()
            assert not button.icon().isNull()
            assert bar.file_action_group.rect().contains(button.geometry())
        edit_menu = next(menu for menu in _top_menus(shell) if menu.title() == "Editar")
        assert edit_menu.actions().count(shell.action_paste) == 1
        assert shell.action_paste not in shell.main_toolbar.actions()
        if width == 1440:
            shell.resize(1440, 900)
            QApplication.processEvents()
            screenshot = tmp_path / "compact-header-light-1440x900.png"
            assert shell.grab().save(str(screenshot))
            assert screenshot.stat().st_size > 0


@pytest.mark.parametrize("theme", ["light", "dark"])
@pytest.mark.parametrize("size", [(1440, 900), (980, 600), (900, 560)])
def test_menu_bar_stays_compact_and_usable_across_themes_and_sizes(
    shell, theme, size
) -> None:
    shell.current_theme = theme
    shell._apply_theme()
    shell.resize(*size)
    QApplication.processEvents()
    menubar = shell.menuBar()
    assert menubar.isVisible()
    assert menubar.height() <= MAX_MENU_BAR_HEIGHT
    assert [menu.title() for menu in _top_menus(shell)] == list(MENU_TITLES)
    for action in menubar.actions():
        geometry = menubar.actionGeometry(action)
        assert geometry.width() > 20 and geometry.height() > 12
    assert shell.app_bar.isVisible()
    assert shell.tool_rail.isVisible()
    assert shell.tabs.width() >= 400


def test_hiding_menu_bar_keeps_hamburger_access_to_the_same_toggle(shell) -> None:
    menus = _top_menus(shell)
    action = shell.action_show_menu_bar
    assert action.text() == "Mostrar barra de herramientas"
    assert shell.view_menu.actions().count(action) == 1
    assert sum(a.text() == action.text() for a in shell.view_menu.actions()) == 1

    action.trigger()
    QApplication.processEvents()
    assert not shell.menuBar().isVisible()
    assert not action.isChecked()
    assert shell.app_bar.isVisible()
    assert shell.main_toolbar.isHidden()
    assert load_ui_preferences(shell._settings).show_menu_bar is False

    # La hamburguesa abre los mismos QMenus; Ver contiene la misma QAction
    # checkable, de modo que no depende de la QMenuBar oculta para restaurarse.
    assert shell.app_bar.menu_button.isVisible()
    shell._open_main_menu_popup()
    QApplication.processEvents()
    aggregated = [
        item.menu()
        for item in shell._menu_aggregator.actions()
        if item.menu() is not None
    ]
    assert aggregated == menus
    assert shell.view_menu in aggregated
    assert action in shell.view_menu.actions()
    assert shell._menu_aggregator.isVisible()
    action.trigger()
    QApplication.processEvents()
    assert shell.menuBar().isVisible()
    assert action.isChecked()
    assert shell.main_toolbar.isHidden()
    assert load_ui_preferences(shell._settings).show_menu_bar is True
    shell._menu_aggregator.close()
    QApplication.processEvents()


def test_new_menu_bar_default_migrates_once_over_legacy_preferences(
    tmp_path, monkeypatch
) -> None:
    """Un INI previo a la beta con el modo compacto persistido migra una vez."""
    legacy = {
        "ui/theme": "dark",
        "ui/shell_defaults_version": 1,
        "ui/show_menu_bar": False,
        "ui/show_main_toolbar": True,
    }
    win = _open_window(tmp_path, monkeypatch, legacy)
    try:
        assert win.current_theme == "dark"
        assert win.menuBar().isVisible()
        assert not win.main_toolbar.isVisible()
        assert win.action_show_menu_bar.isChecked()
        assert win.action_show_menu_bar.text() == "Mostrar barra de herramientas"

        # Una elección explícita posterior se respeta (y sella la migración).
        win.action_show_menu_bar.trigger()
        QApplication.processEvents()
        preferences = load_ui_preferences(win._settings)
        assert preferences.show_menu_bar is False
        assert preferences.show_main_toolbar is False
        assert (
            win._settings.value("ui/shell_defaults_version") == SHELL_DEFAULTS_VERSION
        )
    finally:
        _close_window(win)

    reopened = ChemusonWindow()
    reopened.show()
    QApplication.processEvents()
    try:
        assert not reopened.menuBar().isVisible()
        assert reopened.main_toolbar.isHidden()
        assert not reopened.action_show_menu_bar.isChecked()
    finally:
        _close_window(reopened)


# ---------------------------------------------------------------------------
# A2 — la QToolBar histórica permanece oculta
# ---------------------------------------------------------------------------


def test_main_toolbar_exposes_only_new_open_save_with_existing_actions(shell) -> None:
    toolbar = shell.main_toolbar
    assert toolbar.isHidden()
    assert not toolbar.toggleViewAction().isEnabled()
    actions = [action for action in toolbar.actions() if not action.isSeparator()]
    assert actions == [shell.action_new, shell.action_open, shell.action_save]
    # Son exactamente las QActions de los menús (misma identidad, mismos
    # atajos y mismos controladores): no hay duplicados.
    file_menu = _top_menus(shell)[0]
    assert [a for a in (shell.action_new, shell.action_open, shell.action_save)] == [
        shell.action_new,
        shell.action_open,
        shell.action_save,
    ]
    for action in (shell.action_new, shell.action_open, shell.action_save):
        assert action in file_menu.actions()
    # Las acciones retiradas de la barra siguen en sus rutas actuales.
    for removed in (
        shell.action_undo,
        shell.action_redo,
        shell.action_rotate_left,
        shell.action_rotate_right,
        shell.action_flip_horizontal,
        shell.action_flip_vertical,
        shell.action_clean_2d,
        shell.action_draw_smiles,
    ):
        assert removed not in toolbar.actions()
    edit_menu = _top_menus(shell)[1]
    assert shell.action_undo in edit_menu.actions()
    assert shell.action_redo in edit_menu.actions()
    structure_menu = _top_menus(shell)[3]
    assert shell.action_clean_2d in structure_menu.actions()
    # Undo/Redo: única representación visible en la app bar.
    assert shell.app_bar.undo_button.defaultAction() is shell.action_undo
    assert shell.app_bar.redo_button.defaultAction() is shell.action_redo
    assert shell.app_bar.undo_button.isVisible()
    assert shell.app_bar.redo_button.isVisible()


@pytest.mark.parametrize("theme", ["light", "dark"])
@pytest.mark.parametrize("size", [(1440, 900), (980, 600), (900, 560)])
def test_toolbar_label_toggles_only_menu_bar_across_themes_and_sizes(
    shell, theme, size
) -> None:
    shell.current_theme = theme
    shell._apply_theme()
    shell.resize(*size)
    QApplication.processEvents()
    action = shell.action_show_menu_bar
    toolbar = shell.main_toolbar
    initial_menu_geometry = shell.menuBar().geometry()
    initial_app_bar_geometry = shell.app_bar.geometry()
    assert action.isChecked()
    assert shell.menuBar().isVisible()

    action.trigger()
    QApplication.processEvents()
    assert not shell.menuBar().isVisible()
    assert not action.isChecked()
    assert shell.app_bar.isVisible()
    assert shell.main_toolbar.isHidden()
    assert not toolbar.toggleViewAction().isEnabled()

    action.trigger()
    QApplication.processEvents()
    assert shell.menuBar().isVisible()
    assert action.isChecked()
    assert shell.menuBar().geometry() == initial_menu_geometry
    assert shell.app_bar.geometry() == initial_app_bar_geometry
    assert shell.main_toolbar.isHidden()
    assert shell.tool_rail.isVisible()
    assert shell.side_panel.isHidden()


def test_obsolete_toolbar_preference_never_reveals_legacy_toolbar(
    tmp_path, monkeypatch
) -> None:
    values = {
        "ui/shell_defaults_version": SHELL_DEFAULTS_VERSION,
        "ui/show_menu_bar": True,
        "ui/show_main_toolbar": True,
    }
    _pin_settings(tmp_path, monkeypatch, values)
    win = ChemusonWindow()
    win.show()
    QApplication.processEvents()
    try:
        assert win.menuBar().isVisible()
        assert win.action_show_menu_bar.isChecked()
        assert win.main_toolbar.isHidden()
        assert not win.main_toolbar.toggleViewAction().isEnabled()
        win.action_show_menu_bar.trigger()
        QApplication.processEvents()
        assert not win.menuBar().isVisible()
        assert win.main_toolbar.isHidden()
        assert load_ui_preferences(win._settings).show_main_toolbar is True
    finally:
        _close_window(win)

    # La elección se conserva al reiniciar, mientras el valor legado true de
    # la toolbar no la resucita; la misma QAction sigue alcanzable desde Ver
    # en la hamburguesa y puede restaurar la fila.
    reopened = ChemusonWindow()
    reopened.show()
    QApplication.processEvents()
    try:
        assert not reopened.menuBar().isVisible()
        assert not reopened.action_show_menu_bar.isChecked()
        assert reopened.main_toolbar.isHidden()
        assert reopened.app_bar.menu_button.isVisible()
        reopened._open_main_menu_popup()
        QApplication.processEvents()
        assert reopened.view_menu in [
            item.menu()
            for item in reopened._menu_aggregator.actions()
            if item.menu() is not None
        ]
        assert reopened.action_show_menu_bar in reopened.view_menu.actions()
        reopened.action_show_menu_bar.trigger()
        QApplication.processEvents()
        assert reopened.menuBar().isVisible()
        assert reopened.action_show_menu_bar.isChecked()
        assert reopened.main_toolbar.isHidden()
        assert load_ui_preferences(reopened._settings).show_menu_bar is True
    finally:
        _close_window(reopened)


# ---------------------------------------------------------------------------
# A3 — geometría de pestaña estable (X sin desplazamiento)
# ---------------------------------------------------------------------------


def test_document_tab_bar_keeps_identical_geometry_clean_and_dirty() -> None:
    """Prueba directa: mismo título, sólo cambia dirty → geometría idéntica."""
    bar = DocumentTabBar()
    bar.resize(600, 40)
    bar.show()
    QApplication.processEvents()
    title = "Cadena de ejemplo con nombre largo"
    bar.sync_tabs([title], [False], 0)
    QApplication.processEvents()
    clean = _side_snapshot(bar, 0)
    assert _close_fits_in_tab(clean)

    bar.set_tab(0, title, True)
    QApplication.processEvents()
    assert bar._side_button(0).is_dirty_shown()
    assert _side_snapshot(bar, 0) == clean

    bar.set_tab(0, title, False)
    QApplication.processEvents()
    assert not bar._side_button(0).is_dirty_shown()
    assert _side_snapshot(bar, 0) == clean
    bar.close()
    QApplication.processEvents()


def test_dirty_state_does_not_change_tab_or_close_button_geometry(shell) -> None:
    bar = shell.app_bar.tab_bar
    canvas = shell.canvas
    assert bar.tab_count() == 1
    clean_snapshot = _side_snapshot(bar, 0)
    _, side_size, close_size, close_pos, close_in_tab, right_gap = clean_snapshot
    assert (side_size.width(), side_size.height()) == (_SIDE_WIDTH, _CLOSE_BTN_SIZE)
    assert (close_size.width(), close_size.height()) == (_CLOSE_BTN_SIZE,) * 2
    assert (close_pos.x(), close_pos.y()) == (_SIDE_MARGIN_X + _DOT_SIZE + _SIDE_SPACING, 0)
    assert _close_fits_in_tab(clean_snapshot)
    assert 0 <= right_gap <= 8

    _mark_dirty(canvas, "C")
    assert bar._side_button(0).is_dirty_shown()
    dirty_snapshot = _side_snapshot(bar, 0)
    assert _invariants(dirty_snapshot) == _invariants(clean_snapshot)
    assert _close_fits_in_tab(dirty_snapshot)

    canvas.undo_stack.setClean()
    QApplication.processEvents()
    assert not bar._side_button(0).is_dirty_shown()
    assert _invariants(_side_snapshot(bar, 0)) == _invariants(clean_snapshot)


def test_dirty_dot_paints_only_when_dirty_and_close_icon_always_visible(shell) -> None:
    bar = shell.app_bar.tab_bar
    canvas = shell.canvas
    clean_snapshot = _side_snapshot(bar, 0)
    for theme in ("light", "dark"):
        shell.current_theme = theme
        shell._apply_theme()
        QApplication.processEvents()
        accent = bar._icons.theme_color(theme, "accent")
        icon_tint = bar._icons.theme_color(theme, "icon")
        side = bar._side_button(0)
        side.close_button.clearFocus()
        QApplication.processEvents()

        clean_colors = _rendered_colors(side)
        assert clean_colors[accent] == 0
        assert clean_colors[icon_tint] > 0

        _mark_dirty(canvas, "N")
        dirty_colors = _rendered_colors(side)
        assert dirty_colors[accent] > 0
        # La X sigue pintada con el tinte del tema y en el mismo sitio.
        assert dirty_colors[icon_tint] == clean_colors[icon_tint]
        assert _invariants(_side_snapshot(bar, 0)) == _invariants(clean_snapshot)

        canvas.undo_stack.setClean()
        QApplication.processEvents()
        assert _rendered_colors(side)[accent] == 0


def test_geometry_stable_for_clean_dirty_long_and_reordered_tabs(shell) -> None:
    bar = shell.app_bar.tab_bar
    titles = ("Sin título", "un-titulo-muy-largo-para-forzar-eliremo-de-pestaña")
    shell._on_file_new()
    QApplication.processEvents()
    first, second = (shell._canvas_from_tab_index(i) for i in range(2))
    shell._tab_manager.tab_titles[first] = titles[0]
    shell._tab_manager.tab_titles[second] = titles[1]
    shell._sync_app_bar_tabs()
    QApplication.processEvents()
    _mark_dirty(second, "O")
    snapshots_clean = [_side_snapshot(bar, i) for i in range(2)]

    shell.app_bar.tabMoved.emit(1, 0)
    QApplication.processEvents()
    assert shell.tabs.indexOf(second) == 0
    moved = [_side_snapshot(bar, i) for i in range(2)]
    # El conjunto de geometrías laterales es el mismo tras reordenar: la X no
    # se desplaza ni se sale de su pestaña (lo único que cambia es el ancho
    # del título de la pestaña).
    assert [_invariants(snapshot) for snapshot in moved] == [
        _invariants(snapshot) for snapshot in snapshots_clean
    ]
    for snapshot in moved:
        assert _close_fits_in_tab(snapshot)


def test_close_button_still_targets_its_document_when_dirty(shell, monkeypatch) -> None:
    bar = shell.app_bar.tab_bar
    target = shell.canvas
    _mark_dirty(target, "S")
    shell._on_file_new()
    QApplication.processEvents()
    remaining = shell.canvas
    questions: list[str] = []
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(
            lambda *_args: (questions.append("asked") or QMessageBox.StandardButton.Discard)
        ),
    )
    side = bar._side_button(shell.tabs.indexOf(target))
    assert side is not None and side.close_button.isVisible()
    QTest.mouseClick(side.close_button, Qt.MouseButton.LeftButton)
    QApplication.processEvents()

    assert questions == ["asked"]
    assert _documents(shell) == [remaining]


def test_tab_geometry_survives_theme_and_window_resizes(shell) -> None:
    bar = shell.app_bar.tab_bar
    _mark_dirty(shell.canvas, "F")
    reference = _side_snapshot(bar, 0)
    for theme in ("light", "dark"):
        shell.current_theme = theme
        shell._apply_theme()
        for size in ((1440, 900), (980, 600), (900, 560)):
            shell.resize(*size)
            QApplication.processEvents()
            snapshot = _side_snapshot(bar, 0)
            assert _invariants(snapshot) == _invariants(reference)
            assert _close_fits_in_tab(snapshot)
