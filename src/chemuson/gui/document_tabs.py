"""Pestañas de documento modernas de la barra de aplicación (Fase 3).

``DocumentTabBar`` es un **espejo de solo-lectura** del ``QTabWidget`` de
documentos de :class:`~chemuson.gui.main_window.ChemusonWindow` (fuente
única de verdad: ``CanvasTabManager`` + ``QUndoStack``): no guarda estado
propio de documento; la ventana lo resincroniza con:

- ``sync_tabs(titles, dirty, current)`` — resync completo tras cambios
  estructurales (crear/desechar pestaña);
- ``set_tab(index, title, dirty)`` — actualización puntual de título y
  suciedad (contrato ``update_tab_title``, estado real
  ``undo_stack.isClean()``);
- ``select(index)`` — selección de la pestaña activa.

El componente contiene un ``QTabBar`` real (elide right, scroll buttons,
``movable``) y el botón ``+`` de documento nuevo. Las interacciones (clic,
cierre, ``+``, drag) emiten señales que la ventana conecta al flujo
existente (``setCurrentIndex``, ``_on_tab_close_requested``,
``action_new.trigger()``, ``moveTab``). Lenguaje visual del spike aprobado
(``docs/ui-modernization/pyqt6-spike``): icono de documento, punto de
suciedad ``accent``, cierre discreto, ``+`` dashed; QSS por tokens en
``theme/qss.py`` (``#docTabs``).
"""

from __future__ import annotations

from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtWidgets import QHBoxLayout, QFrame, QTabBar, QToolButton, QWidget

from chemuson.gui.theme.icon_provider import IconProvider

__all__ = [
    "DocumentTabBar",
]

_TAB_ICON_SIZE = 14
_CLOSE_ICON_SIZE = 11
_PLUS_ICON_SIZE = 14
_DOT_SIZE = 7
_CLOSE_BTN_SIZE = 18
_NEW_BTN_SIZE = 28


class _TabSideButton(QWidget):
    """Widget derecho de cada pestaña: punto de suciedad + botón cerrar."""

    def __init__(
        self,
        tint: str,
        icon_provider: IconProvider,
        parent: QWidget | None = None,
    ):
        super().__init__(parent)
        self.setObjectName("docTabSide")
        self._tint = tint
        self._icons = icon_provider
        self._bar: DocumentTabBar | None = None

        layout = QHBoxLayout(self)
        layout.setContentsMargins(2, 0, 2, 0)
        layout.setSpacing(5)

        self.dirty_dot = QWidget(self)
        self.dirty_dot.setObjectName("dirtyDot")
        self.dirty_dot.setFixedSize(_DOT_SIZE, _DOT_SIZE)
        self.dirty_dot.hide()

        self.close_button = QToolButton(self)
        self.close_button.setProperty("tabClose", "true")
        self.close_button.setFixedSize(_CLOSE_BTN_SIZE, _CLOSE_BTN_SIZE)
        self.close_button.setIconSize(QSize(_CLOSE_ICON_SIZE, _CLOSE_ICON_SIZE))
        self.close_button.setToolTip("Cerrar pestaña")
        self.close_button.clicked.connect(lambda _checked=False: self._request_close())

        layout.addWidget(self.dirty_dot)
        layout.addWidget(self.close_button)

    def attach(self, bar: "DocumentTabBar") -> None:
        """Vincula el botón a la barra que lo contiene (para cerrar)."""
        self._bar = bar

    def _request_close(self) -> None:
        if self._bar is not None:
            self._bar._close_tab_for_side(self)

    def set_dirty(self, dirty: bool) -> None:
        self.dirty_dot.setVisible(bool(dirty))

    def is_dirty_shown(self) -> bool:
        # ``isVisibleTo``: el estado de visibilidad propio, independiente de
        # si la cadena de padres está mostrada (tests offscreen sin show).
        return self.dirty_dot.isVisibleTo(self)

    def refresh_icons(self, tint: str) -> None:
        self._tint = tint
        self.close_button.setIcon(self._icons.icon("x", tint, _CLOSE_ICON_SIZE))


class DocumentTabBar(QFrame):
    """Barra de pestañas de documento integrada en la app bar.

    Espejo pasivo del ``QTabWidget`` de la ventana: la ventana decide qué
    pestañas existen (``sync_tabs``) y qué están sucias (``set_tab``); el
    widget solo refleja eso y reemite la interacción del usuario:

    - :signal:`tabActivated(int)` — clic en una pestaña (usuario).
    - :signal:`newDocumentRequested()` — botón ``+``.
    - :signal:`closeRequested(int)` — botón ``x`` de una pestaña.
    - :signal:`tabMoved(int, int)` — drag/reorder.
    """

    tabActivated = pyqtSignal(int)
    newDocumentRequested = pyqtSignal()
    closeRequested = pyqtSignal(int)
    tabMoved = pyqtSignal(int, int)

    def __init__(
        self,
        tint: str | None = None,
        parent: QWidget | None = None,
    ):
        super().__init__(parent)
        self.setObjectName("docTabsWrap")

        self._icons = IconProvider()
        self._tint = tint or self._icons.theme_color("light", "icon")

        self.tab_bar = QTabBar(self)
        self.tab_bar.setObjectName("docTabs")
        self.tab_bar.setExpanding(False)
        self.tab_bar.setDrawBase(False)
        self.tab_bar.setElideMode(Qt.TextElideMode.ElideRight)
        self.tab_bar.setMovable(True)
        self.tab_bar.setUsesScrollButtons(True)
        self.tab_bar.setFocusPolicy(Qt.FocusPolicy.TabFocus)

        self.new_button = QToolButton(self)
        self.new_button.setObjectName("tabNewBtn")
        self.new_button.setFixedSize(_NEW_BTN_SIZE, _NEW_BTN_SIZE)
        self.new_button.setIconSize(QSize(_PLUS_ICON_SIZE, _PLUS_ICON_SIZE))
        self.new_button.setToolTip("Nuevo documento (Ctrl+N)")
        self.new_button.clicked.connect(self.newDocumentRequested.emit)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(self.tab_bar, 1)
        layout.addWidget(self.new_button)

        self.tab_bar.tabBarClicked.connect(self._on_tab_clicked)
        self.tab_bar.tabCloseRequested.connect(self.closeRequested)
        self.tab_bar.tabMoved.connect(self.tabMoved)
        self.refresh_icons()

    # ------------------------------------------------------------------
    # Iconos (theme-aware vía IconProvider: el color va en la caché)
    # ------------------------------------------------------------------
    def refresh_icons(self, tint: str | None = None) -> None:
        """Re-tinta los iconos de pestañas (por defecto, tinte actual)."""
        if tint is not None:
            self._tint = tint
        for index in range(self.tab_bar.count()):
            self.tab_bar.setTabIcon(
                index, self._icons.icon("doc", self._tint, _TAB_ICON_SIZE)
            )
            side = self._side_button(index)
            if side is not None:
                side.refresh_icons(self._tint)
        self.new_button.setIcon(self._icons.icon("plus", self._tint, _PLUS_ICON_SIZE))

    # ------------------------------------------------------------------
    # Sincronización desde el QTabWidget (fuente de verdad)
    # ------------------------------------------------------------------
    def sync_tabs(
        self,
        titles: list[str],
        dirty: list[bool],
        current: int,
    ) -> None:
        """Resincroniza el espejo completo con el estado de la ventana."""
        bar = self.tab_bar
        bar.blockSignals(True)
        try:
            for index in range(bar.count() - 1, -1, -1):
                bar.removeTab(index)
            for position, (title, is_dirty) in enumerate(zip(titles, dirty)):
                self._add_tab(position, title, is_dirty)
            if 0 <= current < bar.count():
                bar.setCurrentIndex(current)
        finally:
            bar.blockSignals(False)

    def set_tab(self, index: int, title: str, dirty: bool) -> None:
        """Actualiza título y punto de suciedad de una pestaña del espejo."""
        if not 0 <= index < self.tab_bar.count():
            return
        self.tab_bar.setTabText(index, title)
        side = self._side_button(index)
        if side is not None:
            side.set_dirty(dirty)

    def select(self, index: int) -> None:
        """Selecciona la pestaña activa del espejo (idempotente)."""
        if 0 <= index < self.tab_bar.count():
            self.tab_bar.setCurrentIndex(index)

    def tab_count(self) -> int:
        return self.tab_bar.count()

    def current_index(self) -> int:
        return self.tab_bar.currentIndex()

    # ------------------------------------------------------------------
    # Internos
    # ------------------------------------------------------------------
    def _add_tab(self, position: int, title: str, dirty: bool) -> int:
        icon = self._icons.icon("doc", self._tint, _TAB_ICON_SIZE)
        index = self.tab_bar.insertTab(position, icon, title)
        side = _TabSideButton(self._tint, self._icons, self)
        side.attach(self)
        self.tab_bar.setTabButton(index, QTabBar.ButtonPosition.RightSide, side)
        side.set_dirty(dirty)
        return index

    def _side_button(self, index: int) -> _TabSideButton | None:
        widget = self.tab_bar.tabButton(index, QTabBar.ButtonPosition.RightSide)
        if isinstance(widget, _TabSideButton):
            return widget
        return None

    def _close_tab_for_side(self, side: _TabSideButton) -> None:
        """Resuelve el índice actual de la pestaña del botón y la cierra.

        El índice se resuelve *en el momento del clic*: el botón conserva
        su asociación con la pestaña a través de reordenamientos y
        remociones de otras pestañas (el índice capturado quedaría
        desfasado, el rectángulo no).
        """
        if not side.isVisible():
            return
        target = side.mapTo(self.tab_bar, side.rect().center())
        for index in range(self.tab_bar.count()):
            if self.tab_bar.tabRect(index).contains(target):
                self.tab_bar.tabCloseRequested.emit(index)
                return

    def _on_tab_clicked(self, index: int) -> None:
        self.tabActivated.emit(index)
