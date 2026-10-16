"""Correcciones de shell de la beta: AppBar de dos franjas y pestañas estables.

OpenSpec: ``fix-beta-tab-close-and-toolbar-visibility`` (Fase A de la
corrección solicitada tras la aceptación visual).

Cubre:
- A1: la barra de menús tradicional (Archivo/Editar/Ver/Estructura/Reacción/
  Ayuda) es visible por defecto y compacta, reutiliza los mismos ``QMenu`` y
  ``QAction``, puede ocultarse desde Ver conservando el acceso hamburguesa, y
  el default nuevo migra una sola vez sobre preferencias previas;
- A2: la barra superior opcional conserva Nuevo/Abrir/Guardar, no aparece por
  defecto ni duplica Pegar; la opción de Ver controla su visibilidad;
- A3: pasar de pestaña limpia a modificada no altera la geometría de la
  pestaña ni la posición de la X; el punto dirty sólo se pinta cuando toca y
  el cierre sigue apuntando al documento correcto.
"""

from __future__ import annotations

from collections import Counter

import pytest
from PyQt6.QtCore import QPoint, QRect, QSettings, Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication, QMessageBox, QToolButton

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
#: Techo de altura de la barra de herramientas de tres botones (baseline: 65 px).
MAX_TOOLBAR_HEIGHT = 40


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


def test_default_header_has_one_compact_action_row_and_shared_file_actions(shell, tmp_path) -> None:
    from PyQt6.QtWidgets import QLabel

    for width in (1440, 980, 900):
        shell.resize(width, 600)
        QApplication.processEvents()
        bar = shell.app_bar
        assert shell.menuBar().isVisible()
        assert bar.height() == 54
        assert shell.main_toolbar.isHidden()
        assert not shell.action_show_main_toolbar.isChecked()
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


def test_hiding_menu_bar_keeps_hamburger_access_to_the_same_menus(shell) -> None:
    menus = _top_menus(shell)
    shell.action_show_menu_bar.trigger()
    QApplication.processEvents()
    assert not shell.menuBar().isVisible()
    assert not shell.action_show_menu_bar.isChecked()
    assert load_ui_preferences(shell._settings).show_menu_bar is False

    # El acceso hamburguesa sigue ahí y abre exactamente los mismos menús.
    assert shell.app_bar.menu_button.isVisible()
    shell._open_main_menu_popup()
    QApplication.processEvents()
    aggregated = [
        action.menu()
        for action in shell._menu_aggregator.actions()
        if action.menu() is not None
    ]
    assert aggregated == menus
    assert shell._menu_aggregator.isVisible()
    shell._menu_aggregator.close()
    QApplication.processEvents()

    shell.action_show_menu_bar.trigger()
    QApplication.processEvents()
    assert shell.menuBar().isVisible()
    assert load_ui_preferences(shell._settings).show_menu_bar is True


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
        assert not win.action_show_main_toolbar.isChecked()

        # Una elección explícita posterior se respeta (y sella la migración).
        win.action_show_menu_bar.trigger()
        win.action_show_main_toolbar.trigger()
        QApplication.processEvents()
        preferences = load_ui_preferences(win._settings)
        assert preferences.show_menu_bar is False
        assert preferences.show_main_toolbar is True
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
        assert reopened.main_toolbar.isVisible()
        assert not reopened.action_show_menu_bar.isChecked()
        assert reopened.action_show_main_toolbar.isChecked()
    finally:
        _close_window(reopened)


# ---------------------------------------------------------------------------
# A2 — barra de herramientas de tres botones
# ---------------------------------------------------------------------------


def test_main_toolbar_exposes_only_new_open_save_with_existing_actions(shell) -> None:
    toolbar = shell.main_toolbar
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
def test_main_toolbar_is_small_uniform_and_theme_stable(shell, theme, size) -> None:
    shell.action_show_main_toolbar.trigger()
    QApplication.processEvents()
    shell.current_theme = theme
    shell._apply_theme()
    shell.resize(*size)
    QApplication.processEvents()
    toolbar = shell.main_toolbar
    assert toolbar.isVisible()
    assert toolbar.height() <= MAX_TOOLBAR_HEIGHT
    buttons = [
        button
        for button in toolbar.findChildren(QToolButton)
        if button.defaultAction() is not None and button.isVisible()
    ]
    assert len(buttons) == 3
    sizes = {button.size() for button in buttons}
    assert len(sizes) == 1, sizes
    button_size = next(iter(sizes))
    assert 20 <= button_size.width() <= 40
    assert 20 <= button_size.height() <= 40
    for button in buttons:
        assert not button.icon().isNull()
        assert button.iconSize().width() == 18
        # Los tres botones caben dentro de la barra (sin recortes).
        assert toolbar.rect().contains(button.geometry())
    # La fila superior sigue siendo compacta en conjunto.
    menubar_bottom = shell.menuBar().mapTo(shell, QPoint(0, shell.menuBar().height())).y()
    assert toolbar.mapTo(shell, QPoint(0, 0)).y() >= menubar_bottom - 1
    assert shell.app_bar.isVisible()


def test_toolbar_toggle_controls_only_the_three_button_bar(shell) -> None:
    toolbar = shell.main_toolbar
    assert not toolbar.isVisible()
    before = (
        shell.menuBar().isVisible(),
        shell.app_bar.isVisible(),
        shell.tool_rail.isVisible(),
        shell.side_panel.isHidden(),
    )
    shell.action_show_main_toolbar.trigger()
    QApplication.processEvents()
    assert toolbar.isVisible()
    assert shell.action_show_main_toolbar.isChecked()
    assert (
        shell.menuBar().isVisible(),
        shell.app_bar.isVisible(),
        shell.tool_rail.isVisible(),
        shell.side_panel.isHidden(),
    ) == before
    assert load_ui_preferences(shell._settings).show_main_toolbar is True

    shell.action_show_main_toolbar.trigger()
    QApplication.processEvents()
    assert not toolbar.isVisible()
    assert [a for a in toolbar.actions() if not a.isSeparator()] == [
        shell.action_new,
        shell.action_open,
        shell.action_save,
    ]
    assert load_ui_preferences(shell._settings).show_main_toolbar is False


def test_toolbar_content_and_visibility_survive_reopening(shell) -> None:
    shell.action_show_main_toolbar.trigger()
    QApplication.processEvents()
    assert shell.main_toolbar.isVisible()
    second = ChemusonWindow()
    second.show()
    QApplication.processEvents()
    try:
        assert second.main_toolbar.isVisible()
        assert [a for a in second.main_toolbar.actions() if not a.isSeparator()] == [
            second.action_new,
            second.action_open,
            second.action_save,
        ]
    finally:
        _close_window(second)


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
