"""Barra de aplicación de Chemuson (Fase 3 del plan de modernización).

``AppBar`` es la franja compacta de documentos y acciones: QActions existentes
para Nuevo/Abrir/Guardar/Pegar, ``DocumentTabBar`` (espejo del ``QTabWidget``
de documentos; sin botón Nuevo duplicado en la ventana real), píldora de
búsqueda (``SearchPill``) que abre la ``CommandPalette`` y controles
undo/redo/tema/ajustes que sostienen las QActions existentes. No presenta
marca, icono de matraz ni versión.

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

_BAR_BUTTON_SIZE = 32
_BAR_ICON_SIZE = 18
_SEARCH_ICON_SIZE = 14
_PILL_WIDTH = 250


def _tooltip_from_action(action: QAction) -> str:
    """Tooltip unificado (Fase 7) para un botón que sostiene una ``QAction``.

    Convención ``Nombre (shortcut)``: si la ``QAction`` expone un shortcut no
    vacío (``QAction.shortcut().toString()``) se incluye; en caso contrario se
    muestra solo el nombre. El ``setDefaultAction`` del botón adopta el
    tooltip de la acción, pero fijar un tooltip explícito (``ownToolTip``)
    preserva la convención sin duplicar handlers ni atajos.
    """
    name = action.text()
    shortcut = action.shortcut().toString()
    if shortcut and "(" not in name:
        return f"{name} ({shortcut})"
    return name


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
    """Píldora de búsqueda de la app bar: entrada a la paleta de comandos.

    **Fase 6**: la píldora es la entrada visible a la ``CommandPalette`` (Ctrl+P).
    Al pulsarla emite :signal:`activated` que la ventana conecta a la misma
    ``QAction`` de apertura que el atajo global ``Ctrl+P`` (un único camino).
    El badge ``Ctrl P`` se muestra (``self.kbd``) y el tooltip deja de indicar
    "(próximamente)". La píldora **no** es un editor permanente: el campo
    editable pertenece a ``CommandPalette``.
    """

    activated = pyqtSignal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("searchPill")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("Buscar o ejecutar… (Ctrl+P)")
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
        # El badge Ctrl P se muestra (Fase 6): la píldora y Ctrl+P abren la
        # misma paleta de comandos. (Ctrl+K volvió a ser Clean2D quick.)
        self.kbd = KbdHint("Ctrl P", self)
        layout.addWidget(self.icon_label)
        layout.addWidget(self.text_label)
        layout.addStretch(1)
        layout.addWidget(self.kbd)

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
    """Franja compacta con documentos y acciones de archivo existentes.

    No crea QAction ni atajos: recibe las QAction originales de la ventana y
    las muestra mediante ``QToolButton``.
    """

    #: Clic en una pestaña del espejo (la ventana activa esa pestaña).
    tabActivated = pyqtSignal(int)
    #: Botón de cierre de una pestaña (la ventana usa el flujo existente).
    closeRequested = pyqtSignal(int)
    #: Reordenamiento por drag en el espejo (la ventana mueve en el
    #: ``QTabWidget`` vía ``moveTab``).
    tabMoved = pyqtSignal(int, int)
    #: Píldora de búsqueda pulsada (Fase 6: la ventana la conecta a la misma
    #: ``QAction`` de apertura que Ctrl+P).
    searchActivated = pyqtSignal()
    #: Botón hamburguesa: la ventana abre los seis menús originales como popup.
    menuRequested = pyqtSignal()

    def __init__(
        self,
        *,
        new_action: QAction,
        open_action: QAction,
        save_action: QAction,
        paste_action: QAction,
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
        # Este botón conserva acceso a los menús originales cuando el
        # ``QMenuBar`` está oculto; las ``QMenu``, ``QAction`` y atajos siguen
        # siendo los mismos tanto en la fila clásica como en el popup.
        self.menu_button = QToolButton(self)
        self.menu_button.setProperty("appBarBtn", "true")
        self.menu_button.setFixedSize(_BAR_BUTTON_SIZE, _BAR_BUTTON_SIZE)
        self.menu_button.setIconSize(QSize(_BAR_ICON_SIZE, _BAR_ICON_SIZE))
        self.menu_button.setToolTip("Menú (Alt)")
        self.menu_button.clicked.connect(self.menuRequested)
        layout.addWidget(self.menu_button)

        # --- Acciones compactas de archivo -------------------------------
        self.file_action_group = QFrame(self)
        self.file_action_group.setObjectName("appBarFileActions")
        file_layout = QHBoxLayout(self.file_action_group)
        file_layout.setContentsMargins(0, 0, 0, 0)
        file_layout.setSpacing(4)
        self.file_action_buttons: dict[str, QToolButton] = {}
        for name, action in (
            ("new", new_action),
            ("open", open_action),
            ("save", save_action),
            ("paste", paste_action),
        ):
            button = QToolButton(self.file_action_group)
            button.setObjectName(f"appBar{name.title()}Button")
            button.setProperty("appBarBtn", "true")
            button.setFixedSize(_BAR_BUTTON_SIZE, _BAR_BUTTON_SIZE)
            button.setIconSize(QSize(_BAR_ICON_SIZE, _BAR_ICON_SIZE))
            button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
            button.setAccessibleName(action.text())
            button.setDefaultAction(action)
            button.setToolTip(_tooltip_from_action(action))
            self.file_action_buttons[name] = button
            file_layout.addWidget(button)
        layout.addWidget(self.file_action_group)

        # --- Pestañas de documento ---------------------------------------
        self.tab_bar = DocumentTabBar(
            self._tint if self._tint else None,
            self,
            show_new_button=False,
        )
        self.tab_bar.tabActivated.connect(self.tabActivated)
        self.tab_bar.closeRequested.connect(self.closeRequested)
        self.tab_bar.tabMoved.connect(self.tabMoved)
        layout.addWidget(self.tab_bar, 1)

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
        # Fase 7: tooltips uniformes ``Nombre (shortcut)``. Se fijan en la
        # ``QAction`` (fuente de verdad): ``setDefaultAction`` los hereda y
        # sobreviven a los re-sync de Qt cuando cambia el estado (enabled) de
        # la acción (un ``setToolTip`` solo en el botón se perdería). La
        # convención se deriva del atajo real de la acción.
        undo_action.setToolTip(_tooltip_from_action(undo_action))
        redo_action.setToolTip(_tooltip_from_action(redo_action))
        theme_action.setToolTip(_tooltip_from_action(theme_action))
        preferences_action.setToolTip(_tooltip_from_action(preferences_action))

        separator = QFrame(self)
        separator.setObjectName("abarSep")
        separator.setFixedHeight(22)

        layout.addWidget(self.undo_button)
        layout.addWidget(self.redo_button)
        layout.addWidget(separator)
        layout.addWidget(self.theme_button)
        layout.addWidget(self.preferences_button)

        # Estado inicial (la ventana lo re-aplica con su tema en
        # ``_apply_theme``; sin esto la pill nacería sin icono).
        self.refresh_icons("light")

    # ------------------------------------------------------------------
    # Iconos (theme-aware; el color es parte de la clave de caché)
    # ------------------------------------------------------------------
    def refresh_icons(self, theme_name: str) -> None:
        """Re-tinta la pill, el botón de tema y las pestañas.

        Los iconos de undo/redo/ajustes son de las QAction (los refresca
        ``refresh_main_toolbar_icons``); aquí se actualizan los propios de
        la barra.
        """
        resolved = "dark" if theme_name == "dark" else "light"
        icon_tint = self._icons.theme_color(resolved, "icon")
        muted_tint = self._icons.theme_color(resolved, "text3")
        self._tint = icon_tint
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
