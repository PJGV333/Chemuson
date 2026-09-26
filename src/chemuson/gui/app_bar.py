"""Barra de aplicación de Chemuson (Fase 3 del plan de modernización).

``AppBar`` es el shell superior nativo (QtWidgets + QSS de tokens + iconos
SVG del ``IconProvider``): marca (``flask`` + nombre + versión),
``DocumentTabBar`` (espejo del ``QTabWidget`` de documentos), botón ``+``,
píldora de búsqueda **placeholder** de la futura command palette (Fase 6;
no registra atajo: Ctrl+K pertenece hoy a ``action_clean_2d_full``; el hint
visual ``Ctrl K`` queda oculto hasta que la Fase 6 implemente la paleta) y
botones undo/redo/tema/ajustes que **sostienen las QAction existentes** de
la ventana (``setDefaultAction`` → icono, tooltip con atajo y estado
habilitado reales; sin duplicar handlers ni shortcuts).

Lenguaje visual del spike aprobado
(``docs/ui-modernization/pyqt6-spike``): altura 54 px (``METRICS``),
superficie ``surface`` + borde inferior, acento cyan en selected/tab dirty.
QSS por tokens en ``theme/qss.py`` (``#app_bar``).
"""

from __future__ import annotations

from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QToolButton,
    QWidget,
)

from chemuson.gui.document_tabs import DocumentTabBar
from chemuson.gui.theme import METRICS
from chemuson.gui.theme.icon_provider import IconProvider

__all__ = [
    "AppBar",
    "SearchPill",
]

_BRAND_ICON_SIZE = 24
_BAR_BUTTON_SIZE = 32
_BAR_ICON_SIZE = 18
_SEARCH_ICON_SIZE = 14
_PILL_WIDTH = 250


class KbdHint(QLabel):
    """Tecla estilo <kbd> (pista visual; no es un atajo registrado)."""

    def __init__(self, text: str, parent: QWidget | None = None):
        super().__init__(text, parent)
        self.setObjectName("kbdK")
        font = self.font()
        font.setPixelSize(10)
        font.setBold(True)
        self.setFont(font)


class SearchPill(QFrame):
    """Píldora de búsqueda de la app bar.

    **Placeholder de la Fase 6** (command palette): visual según el mockup
    (icono + "Buscar o ejecutar…") pero sin atajo ni lógica en esta fase;
    al pulsarla emite :signal:`activated` para que la Fase 6 conecte la
    paleta sin cambiar este widget.

    El badge ``Ctrl K`` queda **oculto** (``self.kbd``) hasta la Fase 6:
    hoy Ctrl+K pertenece a ``action_clean_2d_full`` y no se cambia ese
    atajo; el hint regresa junto con la command palette.
    """

    activated = pyqtSignal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("searchPill")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("Buscar o ejecutar… (próximamente)")
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        # Ancho fijo (250 px, igual que el spike aprobado): la compresión
        # bajo 980 px la asumen las pestañas (elisión), no la píldora.
        self.setFixedWidth(_PILL_WIDTH)

        self._icons = IconProvider()
        self._tint: str | None = None

        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 6, 10, 6)
        layout.setSpacing(8)
        self.icon_label = QLabel(self)
        self.text_label = QLabel("Buscar o ejecutar…", self)
        self.text_label.setObjectName("searchPillTxt")
        # El badge Ctrl K se mantiene construido (lo reactiva la Fase 6)
        # pero oculto: Ctrl+K sigue siendo el de Clean2D full.
        self.kbd = KbdHint("Ctrl K", self)
        self.kbd.hide()
        layout.addWidget(self.icon_label)
        layout.addWidget(self.text_label)
        layout.addStretch(1)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.activated.emit()
        super().mousePressEvent(event)

    def refresh_icons(self, tint: str | None = None) -> None:
        if tint is not None:
            self._tint = tint
        self.icon_label.setPixmap(
            self._icons.pixmap("search", self._tint, _SEARCH_ICON_SIZE)
        )


class AppBar(QFrame):
    """Barra de aplicación: marca, pestañas de documento y controles.

    No crea QAction ni atajos: recibe los existentes de la ventana
    (``undo_action``, ``redo_action``, ``preferences_action``,
    ``theme_action``) y los muestra mediante ``QToolButton``.
    """

    #: Clic en una pestaña del espejo (la ventana activa esa pestaña).
    tabActivated = pyqtSignal(int)
    #: Botón ``+`` (la ventana dispara la QAction "nuevo" existente).
    newDocumentRequested = pyqtSignal()
    #: Botón de cierre de una pestaña (la ventana usa el flujo existente).
    closeRequested = pyqtSignal(int)
    #: Reordenamiento por drag en el espejo (la ventana mueve en el
    #: ``QTabWidget`` vía ``moveTab``).
    tabMoved = pyqtSignal(int, int)
    #: Píldora de búsqueda pulsada (placeholder; lo conectará la Fase 6).
    searchActivated = pyqtSignal()
    #: Botón hamburguesa: la ventana abre el ``QMenuBar`` (oculto) como popup.
    menuRequested = pyqtSignal()

    def __init__(
        self,
        *,
        version: str,
        undo_action: QAction,
        redo_action: QAction,
        preferences_action: QAction,
        theme_action: QAction,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("app_bar")
        self.setFixedHeight(METRICS["appbarH"])

        self._icons = IconProvider()
        self._tint: str | None = None

        layout = QHBoxLayout(self)
        layout.setContentsMargins(12, 0, 12, 0)
        layout.setSpacing(12)

        # --- Menú (hamburguesa) ----------------------------------------
        # El ``QMenuBar`` histórico queda oculto; este botón es el acceso
        # visible (clic o tecla Alt → ``QMenuBar.popup``). Las ``QMenu``,
        # ``QAction`` y atajos siguen siendo las originales.
        self.menu_button = QToolButton(self)
        self.menu_button.setProperty("appBarBtn", "true")
        self.menu_button.setFixedSize(_BAR_BUTTON_SIZE, _BAR_BUTTON_SIZE)
        self.menu_button.setIconSize(QSize(_BAR_ICON_SIZE, _BAR_ICON_SIZE))
        self.menu_button.setToolTip("Menú (Alt)")
        self.menu_button.clicked.connect(self.menuRequested)
        layout.addWidget(self.menu_button)

        # --- Marca -------------------------------------------------------
        self.brand_label = QLabel(self)
        layout.addWidget(self.brand_label)
        self.brand_name = QLabel("Chemuson", self)
        self.brand_name.setObjectName("appBrandName")
        layout.addWidget(self.brand_name)
        self.version_label = QLabel(version or "", self)
        self.version_label.setObjectName("appVersionPill")
        layout.addWidget(self.version_label)

        # --- Pestañas de documento + botón nuevo -------------------------
        self.tab_bar = DocumentTabBar(self._tint if self._tint else None, self)
        self.tab_bar.tabActivated.connect(self.tabActivated)
        self.tab_bar.newDocumentRequested.connect(self.newDocumentRequested)
        self.tab_bar.closeRequested.connect(self.closeRequested)
        self.tab_bar.tabMoved.connect(self.tabMoved)
        layout.addWidget(self.tab_bar)

        # --- Controles derechos ------------------------------------------
        self.search_pill = SearchPill(self)
        self.search_pill.activated.connect(self.searchActivated)
        layout.addStretch(1)
        layout.addWidget(self.search_pill)

        self.undo_button = QToolButton(self)
        self.redo_button = QToolButton(self)
        self.theme_button = QToolButton(self)
        self.preferences_button = QToolButton(self)
        for button in (
            self.undo_button,
            self.redo_button,
            self.theme_button,
            self.preferences_button,
        ):
            button.setProperty("appBarBtn", "true")
            button.setFixedSize(_BAR_BUTTON_SIZE, _BAR_BUTTON_SIZE)
            button.setIconSize(QSize(_BAR_ICON_SIZE, _BAR_ICON_SIZE))
        self.undo_button.setDefaultAction(undo_action)
        self.redo_button.setDefaultAction(redo_action)
        self.theme_button.setDefaultAction(theme_action)
        self.preferences_button.setDefaultAction(preferences_action)

        separator = QFrame(self)
        separator.setObjectName("abarSep")
        separator.setFixedHeight(22)

        layout.addWidget(self.undo_button)
        layout.addWidget(self.redo_button)
        layout.addWidget(separator)
        layout.addWidget(self.theme_button)
        layout.addWidget(self.preferences_button)

        # Estado inicial (la ventana lo re-aplica con su tema en
        # ``_apply_theme``; sin esto la marca/pill nacerían sin icono).
        self.refresh_icons("light")

    # ------------------------------------------------------------------
    # Iconos (theme-aware; el color es parte de la clave de caché)
    # ------------------------------------------------------------------
    def refresh_icons(self, theme_name: str) -> None:
        """Re-tinta la marca, la pill, el botón de tema y las pestañas.

        Los iconos de undo/redo/ajustes son de las QAction (los refresca
        ``refresh_main_toolbar_icons``); aquí se actualizan los propios de
        la barra.
        """
        resolved = "dark" if theme_name == "dark" else "light"
        accent = self._icons.theme_color(resolved, "accent")
        icon_tint = self._icons.theme_color(resolved, "icon")
        muted_tint = self._icons.theme_color(resolved, "text3")
        self._tint = icon_tint
        self.brand_label.setPixmap(
            self._icons.pixmap("flask", accent, _BRAND_ICON_SIZE)
        )
        self.menu_button.setIcon(
            self._icons.icon("menu", icon_tint, _BAR_ICON_SIZE)
        )
        self.theme_button.setIcon(
            self._icons.icon("sun" if resolved == "dark" else "moon", icon_tint, _BAR_ICON_SIZE)
        )
        self.search_pill.refresh_icons(muted_tint)
        self.tab_bar.refresh_icons(icon_tint)

    # ------------------------------------------------------------------
    # Sincronización de pestañas (delegado al espejo)
    # ------------------------------------------------------------------
    def sync_tabs(
        self,
        titles: list[str],
        dirty: list[bool],
        current: int,
    ) -> None:
        self.tab_bar.sync_tabs(titles, dirty, current)

    def set_tab(self, index: int, title: str, dirty: bool) -> None:
        self.tab_bar.set_tab(index, title, dirty)

    def select(self, index: int) -> None:
        self.tab_bar.select(index)

    def tab_count(self) -> int:
        return self.tab_bar.tab_count()
