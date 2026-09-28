from __future__ import annotations

import pytest
from PyQt6.QtCore import QSettings, QStandardPaths, QPoint, Qt
from PyQt6.QtGui import QFontMetrics
from PyQt6.QtWidgets import QApplication, QDockWidget, QLabel, QMenu

from chemuson.gui.main_window import ChemusonWindow
from chemuson.gui.side_panel import SidePanel
from chemuson.platform.settings import application_settings


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    config_location = QStandardPaths.writableLocation(
        QStandardPaths.StandardLocation.ConfigLocation
    )
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    QSettings.setPath(
        QSettings.Format.NativeFormat,
        QSettings.Scope.UserScope,
        str(tmp_path),
    )
    try:
        yield
    finally:
        QSettings.setPath(
            QSettings.Format.NativeFormat,
            QSettings.Scope.UserScope,
            config_location,
        )


@pytest.fixture
def window():
    instance = ChemusonWindow()
    try:
        yield instance
    finally:
        instance.close()
        QApplication.processEvents()


def _view_menu(window: ChemusonWindow) -> QMenu:
    return next(
        action.menu()
        for action in window.menuBar().actions()
        if action.menu() is not None
        and action.menu().title().replace("&", "").strip().lower() == "ver"
    )


def _dock_pages(window: ChemusonWindow) -> dict[str, QDockWidget]:
    return {
        "inspector": window.inspector_dock,
        "validation": window.validation_dock,
        "properties": window.chemical_properties_dock,
        "templates": window.templates_dock,
        "appearance": window.appearance_dock,
        "spectroscopy": window.spectroscopy_dock,
        "compchem": window.compchem_dock,
    }


def test_side_panel_hosts_the_original_docks_without_duplicate_pages(window) -> None:
    panel = window.side_panel
    docks = _dock_pages(window)

    assert isinstance(panel, SidePanel)
    assert 336 <= panel.width() <= 340
    assert panel.stack.count() == len(docks) == 7
    assert len(window.findChildren(QDockWidget)) == 7

    for key, dock in docks.items():
        assert panel.page_widget(key) is dock
        assert panel.stack.indexOf(dock) >= 0
        assert dock.widget() is panel.page_widget(key).widget()
        assert window.dockWidgetArea(dock) == Qt.DockWidgetArea.NoDockWidgetArea
        assert dock.features() == QDockWidget.DockWidgetFeature.NoDockWidgetFeatures
        assert not dock.isFloating()


def test_side_panel_defaults_to_inspector_and_overflow_selects_in_place(window) -> None:
    panel = window.side_panel
    assert panel.active_page_key == "inspector"
    assert panel.stack.currentWidget() is window.inspector_dock
    assert tuple(panel.main_tab_buttons) == SidePanel.MAIN_PAGE_KEYS
    assert tuple(action.text() for action in panel.overflow_menu.actions()) == (
        "Espectroscopía",
        "CompChem",
    )

    for key in ("spectroscopy", "compchem"):
        panel.overflow_actions[key].trigger()
        assert panel.active_page_key == key
        assert panel.stack.currentWidget() is _dock_pages(window)[key]
        assert panel.overflow_button.isChecked()
        assert not panel.isHidden()


def test_validation_controls_fit_within_side_panel_page(window) -> None:
    window.resize(1440, 900)
    window.show()
    window.side_panel.show_page("validation")
    dock = window.validation_dock
    controls = (
        dock.btn_refresh,
        dock.btn_previous,
        dock.btn_next,
        dock.btn_copy_report,
        dock.btn_export_report,
        dock.action_combo,
        dock.btn_apply_correction,
    )
    dock.action_combo.addItem("Corrección")
    dock.btn_apply_correction.setEnabled(True)
    for control in controls:
        control.show()
    QApplication.processEvents()

    panel_rect = window.side_panel.rect()
    for control in controls:
        top_left = control.mapTo(window.side_panel, control.rect().topLeft())
        right = top_left.x() + control.width()
        bottom = top_left.y() + control.height()
        label = (
            control.objectName()
            or getattr(control, "text", lambda: "")()
            or type(control).__name__
        )
        assert (
            panel_rect.left() <= top_left.x()
            and right <= panel_rect.right() + 1
            and panel_rect.top() <= top_left.y()
            and bottom <= panel_rect.bottom() + 1
        ), (
            f"validation control {label} "
            f"extends outside the side panel: {top_left}, {control.size()}"
        )


def test_view_menu_actions_show_and_select_pages_without_dock_toggles(window) -> None:
    view_menu = _view_menu(window)
    actions = {action.text(): action for action in view_menu.actions()}
    labels_and_keys = (
        ("Inspector", "inspector"),
        ("Validación", "validation"),
        ("Propiedades químicas", "properties"),
        ("Plantillas", "templates"),
        ("Apariencia", "appearance"),
        ("Espectroscopía", "spectroscopy"),
        ("CompChem", "compchem"),
    )
    docks = _dock_pages(window)

    for label, key in labels_and_keys:
        assert label in actions
        assert all(actions[label] is not dock.toggleViewAction() for dock in docks.values())
        window.side_panel.set_panel_visible(False)
        actions[label].trigger()
        assert not window.side_panel.isHidden()
        assert window.side_panel.active_page_key == key
        assert window.side_panel.stack.currentWidget() is docks[key]
        assert all(not dock.isFloating() for dock in docks.values())


def test_side_panel_persists_active_page_and_visibility(window) -> None:
    panel = window.side_panel
    panel.show_page("appearance")
    panel.set_panel_visible(False)

    settings = window._settings
    assert settings.value("ui/side_panel/active_tab") == "appearance"
    assert settings.value("ui/side_panel/visible") is False

    window.close()
    restored = ChemusonWindow()
    try:
        assert restored.side_panel.active_page_key == "appearance"
        assert restored.side_panel.isHidden()
    finally:
        restored.close()


def test_invalid_side_panel_settings_fall_back_to_visible_inspector() -> None:
    settings = application_settings()
    settings.setValue("ui/side_panel/active_tab", "unknown")
    settings.setValue("ui/side_panel/visible", "")

    window = ChemusonWindow()
    try:
        assert window.side_panel.active_page_key == "inspector"
        assert not window.side_panel.isHidden()
    finally:
        window.close()


def test_app_bar_menu_reaches_appearance_and_preferences_action_is_preserved(window) -> None:
    assert window.app_bar.preferences_button.defaultAction() is window.action_preferences

    actions = {action.text(): action for action in _view_menu(window).actions()}
    actions["Apariencia"].trigger()

    assert window.side_panel.active_page_key == "appearance"
    assert not window.side_panel.isHidden()


def test_status_bar_keeps_existing_indicators_and_show_message(window) -> None:
    window.resize(1440, 900)
    window.show()
    QApplication.processEvents()
    status_bar = window.statusBar()
    assert status_bar.height() == 34
    assert window._iupac_name_label.x() >= 300
    assert window._total_charge_label.x() > window._iupac_name_label.x()

    window._update_iupac_name_indicator()
    window._update_total_charge_indicator()
    status_bar.showMessage("Mensaje temporal")
    QApplication.processEvents()

    assert status_bar.currentMessage() == "Mensaje temporal"
    assert window._iupac_name_label.text().startswith("Nombre IUPAC:")
    assert window._total_charge_label.text()
    assert window.findChild(QLabel, "statusFormulaLabel") is None


def test_primary_tabs_have_complete_labels_padding_and_separation(window) -> None:
    window.show()
    panel = window.side_panel
    row = panel.tab_row
    viewport = row.scroll_area.viewport()
    minimum_padding_x = 3
    minimum_gap = 3

    for size in ((1440, 900), (980, 600)):
        window.resize(*size)
        QApplication.processEvents()
        previous_right = None
        for key, button in row.main_tab_buttons.items():
            position = button.mapTo(viewport, QPoint(0, 0))
            text_width = QFontMetrics(button.font()).horizontalAdvance(button.text())
            right = position.x() + button.width()

            assert button.isVisible(), f"primary tab {key} must remain visible at {size}"
            assert button.font().pixelSize() >= 10, "tab text must not be reduced further"
            assert button.width() >= text_width + 2 * minimum_padding_x, (
                f"primary tab {key} needs at least {minimum_padding_x}px lateral "
                f"padding around its {text_width}px label"
            )
            assert position.x() >= 0 and right <= viewport.width(), (
                f"primary tab {key} is clipped at {size}: "
                f"x={position.x()}, width={button.width()}, viewport={viewport.width()}"
            )
            if previous_right is not None:
                gap = position.x() - previous_right
                assert gap >= minimum_gap, (
                    f"primary tab {key} needs a visible {minimum_gap}px gap; got {gap}px"
                )
            previous_right = right

        assert 336 <= panel.width() <= 340, (
            f"SidePanel width must stay within 336–340px for readable tabs; "
            f"got {panel.width()}px"
        )


def test_side_panel_remains_usable_in_dark_theme_at_compact_window_size(window) -> None:
    active_page = window.side_panel.active_page_key
    window.resize(980, 600)
    window.show()
    window.toggle_theme(True)
    QApplication.processEvents()

    assert 336 <= window.side_panel.width() <= 340
    assert 300 <= window.side_panel.width() <= 340
    assert window.tabs.width() > 0
    assert window.side_panel.active_page_key == active_page

    window.toggle_theme(False)
    QApplication.processEvents()
    assert 336 <= window.side_panel.width() <= 340
    assert window.statusBar().height() == 34
    assert window.side_panel.active_page_key == active_page
