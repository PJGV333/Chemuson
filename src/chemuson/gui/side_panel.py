"""Modern right-side panel that reuses the historical dock widgets.

The panel owns navigation and presentation only. The existing QDockWidget
instances remain the sole owners of their content, state, signals and actions.
"""

from __future__ import annotations

from collections.abc import Mapping

from PyQt6.QtCore import QPoint, Qt, pyqtSignal
from PyQt6.QtGui import QAction, QFontMetrics
from PyQt6.QtWidgets import (
    QDockWidget,
    QFrame,
    QHBoxLayout,
    QMenu,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemuson.gui.theme import METRICS
from chemuson.platform.settings import (
    SIDE_PANEL_TAB_KEYS,
    SettingsStore,
    SidePanelPreferences,
    application_settings,
    load_side_panel_preferences,
    save_side_panel_preferences,
)


class SideTabRow(QFrame):
    """Scrollable primary tabs with a fixed overflow control."""

    pageRequested = pyqtSignal(str)
    overflowRequested = pyqtSignal()

    MAIN_PAGE_KEYS = (
        "inspector",
        "validation",
        "properties",
        "templates",
        "appearance",
    )
    PAGE_LABELS = {
        "inspector": "Inspector",
        "validation": "Validación",
        "properties": "Propiedades",
        "templates": "Plantillas",
        "appearance": "Apariencia",
        "spectroscopy": "Espectroscopía",
        "compchem": "CompChem",
    }
    OVERFLOW_PAGE_KEYS = ("spectroscopy", "compchem")

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("sideTabRow")
        self.setFixedHeight(METRICS["sideTabH"])
        self.main_tab_buttons: dict[str, QToolButton] = {}

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.scroll_area = QScrollArea(self)
        self.scroll_area.setObjectName("sideTabsScroll")
        self.scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll_area.setWidgetResizable(False)
        self.scroll_area.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.scroll_area.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.scroll_area.setFixedHeight(40)
        self._strip = QWidget()
        self._strip.setObjectName("sideTabsStrip")
        self._strip.setFixedHeight(40)
        strip_layout = QHBoxLayout(self._strip)
        strip_layout.setContentsMargins(2, 0, 2, 0)
        strip_layout.setSpacing(0)

        for key in self.MAIN_PAGE_KEYS:
            button = QToolButton(self._strip)
            button.setObjectName(f"sideTab_{key}")
            button.setProperty("sideTab", True)
            button.setText(self.PAGE_LABELS[key])
            button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
            button.setCheckable(True)
            button.setAutoRaise(True)
            button.setCursor(Qt.CursorShape.PointingHandCursor)
            button.setFixedHeight(40)
            tab_font = button.font()
            tab_font.setPixelSize(METRICS["sideTabFont"])
            button.setFont(tab_font)
            tab_text_width = QFontMetrics(button.font()).horizontalAdvance(button.text())
            button.setFixedWidth(tab_text_width + 4)
            button.clicked.connect(
                lambda _checked=False, page=key: self.pageRequested.emit(page)
            )
            strip_layout.addWidget(button)
            self.main_tab_buttons[key] = button

        self._strip.adjustSize()
        self._strip.setMinimumWidth(strip_layout.sizeHint().width())
        self.scroll_area.setWidget(self._strip)
        layout.addWidget(self.scroll_area, 1)

        self.overflow_button = QToolButton(self)
        self.overflow_button.setObjectName("sidePanelOverflow")
        self.overflow_button.setProperty("sideTab", True)
        self.overflow_button.setText("…")
        self.overflow_button.setToolTip("Más paneles")
        self.overflow_button.setAccessibleName("Más paneles")
        self.overflow_button.setCheckable(True)
        self.overflow_button.setAutoRaise(True)
        self.overflow_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.overflow_button.setFixedSize(36, 40)
        self.overflow_button.clicked.connect(self.overflowRequested)
        layout.addWidget(self.overflow_button)

    def set_active(self, page_key: str) -> None:
        """Reflect a primary or overflow page without storing page state."""
        is_overflow = page_key in self.OVERFLOW_PAGE_KEYS
        for key, button in self.main_tab_buttons.items():
            active = key == page_key
            button.setChecked(active)
            button.setProperty("active", active)
            button.style().unpolish(button)
            button.style().polish(button)

        self.overflow_button.setChecked(is_overflow)
        self.overflow_button.setProperty("active", is_overflow)
        self.overflow_button.style().unpolish(self.overflow_button)
        self.overflow_button.style().polish(self.overflow_button)

        button = self.main_tab_buttons.get(page_key)
        if button is not None:
            self.scroll_area.ensureWidgetVisible(button, 8, 0)


class SidePanel(QFrame):
    """Tabbed presentation surface for the seven existing dock instances."""

    activePageChanged = pyqtSignal(str)
    panelVisibilityChanged = pyqtSignal(bool)

    MAIN_PAGE_KEYS = SideTabRow.MAIN_PAGE_KEYS
    OVERFLOW_PAGE_KEYS = SideTabRow.OVERFLOW_PAGE_KEYS
    PAGE_LABELS = SideTabRow.PAGE_LABELS
    PAGE_KEYS = SIDE_PANEL_TAB_KEYS

    def __init__(
        self,
        docks: Mapping[str, QDockWidget],
        *,
        settings: SettingsStore | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        missing = set(self.PAGE_KEYS) - set(docks)
        extra = set(docks) - set(self.PAGE_KEYS)
        if missing or extra:
            raise ValueError(
                f"SidePanel dock keys mismatch (missing={sorted(missing)}, extra={sorted(extra)})"
            )

        self.setObjectName("sidePanel")
        self.setFixedWidth(METRICS["sideW"])
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self._settings = settings or application_settings()
        self._pages = dict(docks)
        self._overflow_actions: dict[str, QAction] = {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.tab_row = SideTabRow(self)
        self.tab_row.pageRequested.connect(self.show_page)
        self.tab_row.overflowRequested.connect(self._show_overflow_menu)
        layout.addWidget(self.tab_row)

        self.stack = QStackedWidget(self)
        self.stack.setObjectName("sideBody")
        layout.addWidget(self.stack, 1)

        for key in self.PAGE_KEYS:
            dock = self._pages[key]
            dock.setFeatures(QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
            dock.setProperty("sidePanelPage", True)
            title_bar = QWidget(dock)
            title_bar.setFixedHeight(0)
            dock.setTitleBarWidget(title_bar)
            self.stack.addWidget(dock)

        self.overflow_menu = QMenu(self)
        self.overflow_menu.setObjectName("sidePanelOverflowMenu")
        self.overflow_menu.addAction(self._make_overflow_action("spectroscopy"))
        self.overflow_menu.addAction(self._make_overflow_action("compchem"))

        preferences = load_side_panel_preferences(self._settings)
        self.active_page_key = preferences.active_tab
        self.stack.setCurrentWidget(self._pages[self.active_page_key])
        self.tab_row.set_active(self.active_page_key)
        self._sync_overflow_actions()
        self.setVisible(preferences.visible)

    @property
    def main_tab_buttons(self) -> dict[str, QToolButton]:
        """Checkable buttons for the five primary pages."""
        return self.tab_row.main_tab_buttons

    @property
    def overflow_button(self) -> QToolButton:
        return self.tab_row.overflow_button

    @property
    def overflow_actions(self) -> dict[str, QAction]:
        """Actions for the two pages reached through the overflow menu."""
        return self._overflow_actions

    def page_widget(self, page_key: str) -> QDockWidget:
        """Return the original dock instance hosted for ``page_key``."""
        try:
            return self._pages[page_key]
        except KeyError as error:
            raise KeyError(f"Unknown side-panel page: {page_key}") from error

    def show_page(self, page_key: str) -> None:
        """Select a page, show the panel and persist the resulting state."""
        if page_key not in self._pages:
            raise KeyError(f"Unknown side-panel page: {page_key}")
        changed = page_key != self.active_page_key
        self.active_page_key = page_key
        self.stack.setCurrentWidget(self._pages[page_key])
        self.tab_row.set_active(page_key)
        self._sync_overflow_actions()
        self.set_panel_visible(True)
        if changed:
            self.activePageChanged.emit(page_key)

    def set_panel_visible(self, visible: bool) -> None:
        """Set panel visibility and persist it without changing the active page."""
        visible = bool(visible)
        was_hidden = self.isHidden()
        self.setVisible(visible)
        save_side_panel_preferences(
            self._settings,
            SidePanelPreferences(active_tab=self.active_page_key, visible=visible),
        )
        if was_hidden != self.isHidden():
            self.panelVisibilityChanged.emit(visible)

    def _make_overflow_action(self, page_key: str) -> QAction:
        action = QAction(self.PAGE_LABELS[page_key], self)
        action.setCheckable(True)
        action.triggered.connect(
            lambda _checked=False, key=page_key: self.show_page(key)
        )
        self._overflow_actions[page_key] = action
        return action

    def _sync_overflow_actions(self) -> None:
        for key, action in self._overflow_actions.items():
            action.setChecked(key == self.active_page_key)

    def _show_overflow_menu(self) -> None:
        point = self.tab_row.overflow_button.mapToGlobal(
            QPoint(0, self.tab_row.overflow_button.height())
        )
        self.overflow_menu.popup(point)
