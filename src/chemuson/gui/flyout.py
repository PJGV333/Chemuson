"""
Flyout reutilizable de paleta para el rail de herramientas (Fase 4).

Componente visual (QFrame frameless, 244 px) que reemplaza visualmente los
``QMenu`` de paleta históricos: cabecera (título + tecla ``Esc``),
cuadrícula de celdas (icono 22 px + etiqueta con word-wrap de 1–3 líneas)
y pie opcional (separador + etiqueta + hasta 3 botones de texto). Se
muestra junto a un botón del rail y se cierra con ``Esc``, clic fuera o
selección de celda.

El flyout **no posee lógica de herramienta**: ejecuta los callbacks que le
pasen (los del toolbar original, ver ``tool_rail.py``) y emite ``closed``
cuando se oculta.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from PyQt6.QtCore import QObject, QEvent, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFontMetrics, QMouseEvent
from PyQt6.QtWidgets import (
    QFrame,
    QGraphicsDropShadowEffect,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemuson.gui.theme.tokens import METRICS, get_tokens

__all__ = [
    "Flyout",
    "FlyoutCell",
    "FlyoutFooter",
    "FlyoutFooterButton",
    "FlyoutItem",
]


@dataclass(frozen=True)
class FlyoutItem:
    """Descriptor de una celda del flyout.

    Args:
        item_id: Identificador estable (para ``on_select``/``active_id``).
        icon: :class:`QIcon` (puede estar vacío: solo etiqueta).
        label: Etiqueta (con word-wrap de 1–3 líneas).
        enabled: Si la celda es interactiva.
        tooltip: Tooltip opcional (por defecto la etiqueta).
        trigger: Callback ejecutado al pulsar la celda (delegación 1:1 al
            owner; p. ej. ``QToolButton.click`` del menú original).
    """

    item_id: str
    icon: object
    label: str
    enabled: bool = True
    tooltip: str = ""
    trigger: Callable[[], None] | None = None


@dataclass(frozen=True)
class FlyoutFooterButton:
    """Botón de texto del pie del flyout (1–3 por flyout).

    Args:
        text: Texto del botón.
        callback: Callback al pulsar (delegación al owner).
        close_after: Cerrar el flyout tras el callback (por defecto sí).
    """

    text: str
    callback: Callable[[], None]
    close_after: bool = True


@dataclass
class FlyoutFooter:
    """Pie opcional del flyout (separador + etiqueta + botones)."""

    text: str = ""
    buttons: list[FlyoutFooterButton] = field(default_factory=list)


def _wrap_text(text: str, fm: QFontMetrics, width: int) -> str:
    """Envuelve ``text`` a ``width`` px (máx. 3 líneas, 2º línea con "...")."""
    words = str(text).split()
    if not words:
        return str(text)
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if fm.horizontalAdvance(candidate) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > 3:
        # Truncar a 3 líneas con el resto indicado.
        last = lines[2]
        while last and fm.horizontalAdvance(f"{last}...") > width:
            last = last.rsplit(" ", 1)[0]
        lines = lines[:2] + [f"{last}..."]
    return "\n".join(lines)


class FlyoutCell(QFrame):
    """Celda de flyout: icono 22 px + etiqueta con word-wrap (1–3 líneas).

    Equivalente al ``FlyoutCell`` del spike PyQt6 (``.flyout-item`` del
    mockup): ``QToolButton`` no envuelve el texto (elide), por esto es un
    ``QFrame`` custom con ``QLabel`` word-wrap.
    """

    clicked = pyqtSignal()

    def __init__(self, icon: object, label: str, column_width: int, parent=None):
        super().__init__(parent)
        self.setProperty("cls", "flyItem")
        tooltip = label
        self.setToolTip(tooltip)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedWidth(column_width)
        # Altura mínima de una línea: sin ella, cuando QGridLayout recalcula
        # mal tras un repopulate, las filas se colapsan a 0 px.
        self.setMinimumHeight(56)
        v = QVBoxLayout(self)
        v.setContentsMargins(2, 8, 2, 7)
        v.setSpacing(5)
        self._icon_label = QLabel()
        self._icon_label.setFixedSize(22, 22)
        self._icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        if icon is not None:
            self._icon_label.setPixmap(icon.pixmap(22, 22))
        v.addWidget(self._icon_label, 0, Qt.AlignmentFlag.AlignHCenter)
        self._label = QLabel(label)
        self._label.setObjectName("flyLbl")
        self._label.setWordWrap(True)
        self._label.setAlignment(
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop
        )
        font = self._label.font()
        font.setPixelSize(10)
        self._label.setFont(font)
        self._label.setText(_wrap_text(label, QFontMetrics(font), column_width - 6))
        v.addWidget(self._label)
        if icon is None:
            self._label.setText(_wrap_text(label, QFontMetrics(font), column_width - 6))

    def set_active(self, on: bool) -> None:
        """Marca/desmarca la celda activa (estado ``active`` QSS)."""
        self.setProperty("active", on)
        style = self.style()
        style.unpolish(self)
        style.polish(self)

    def mousePressEvent(self, e: QMouseEvent) -> None:  # noqa: N802
        if e.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(e)


class Flyout(QFrame):
    """Flyout de paleta: cabecera + cuadrícula + pie opcional.

    Ancho fijo :data:`~chemuson.gui.theme.tokens.METRICS` ``flyoutW`` (244
    px, HTML: ``.flyout { width: 244px; padding: 11px }``). Se muestra con
    :meth:`show_near`; se cierra con ``Esc``, clic fuera o al elegir una
    celda (``close_on_select``). Es hijo de la ventana (no modal).
    """

    closed = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("flyout")
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setFixedWidth(METRICS["flyoutW"])
        self._on_select: Callable[[str], None] | None = None
        self._close_on_select = True
        self._cells: list[FlyoutCell] = []
        self._cell_ids: list[str] = []
        self._foot_buttons: list[QToolButton] = []
        self._close_filter_target: QObject | None = None
        self._shadow = QGraphicsDropShadowEffect(self)
        self._shadow.setBlurRadius(40)
        self._shadow.setOffset(0, 14)
        self._shadow.setColor(QColor(get_tokens("light")["shadow2"]))
        self.setGraphicsEffect(self._shadow)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(11, 11, 11, 11)
        lay.setSpacing(9)

        head = QHBoxLayout()
        head.setSpacing(8)
        self._title = QLabel()
        self._title.setObjectName("flyoutTitle")
        kbd = QLabel("Esc")
        kbd.setObjectName("flyKbd")
        head.addWidget(self._title)
        head.addStretch(1)
        head.addWidget(kbd)
        lay.addLayout(head)

        self._grid = QGridLayout()
        self._grid.setSpacing(6)
        lay.addLayout(self._grid)

        self._foot_sep = QFrame()
        self._foot_sep.setObjectName("flyoutSep")
        self._foot_sep.setFixedHeight(1)
        self._foot_txt = QLabel()
        self._foot_txt.setObjectName("flyFootTxt")
        self._foot_txt.setVisible(False)
        self._foot_row = QHBoxLayout()
        self._foot_row.setSpacing(8)
        self._foot_row.addWidget(self._foot_txt)
        self._foot_row.addStretch(1)
        lay.addWidget(self._foot_sep)
        lay.addLayout(self._foot_row)
        self.hide()

    # ------------------------------------------------------------------
    # API pública
    # ------------------------------------------------------------------
    def set_theme(self, theme_name: str) -> None:
        """Actualiza el color de sombra con los tokens del tema."""
        self._shadow.setColor(QColor(get_tokens(theme_name)["shadow2"]))

    def populate(
        self,
        title: str,
        items: list[FlyoutItem],
        columns: int = 4,
        active_id: str | None = None,
        on_select: Callable[[str], None] | None = None,
        close_on_select: bool = True,
        footer: FlyoutFooter | None = None,
    ) -> None:
        """Reconstruye el flyout (título, celdas y pie).

        Args:
            title: Título de la cabecera.
            items: Celdas a mostrar.
            columns: Columnas de la cuadrícula.
            active_id: Celda a marcar como activa (opcional).
            on_select: Callback con el ``item_id`` elegido.
            close_on_select: Cerrar tras la selección (por defecto sí).
            footer: Pie opcional (etiqueta + 1–3 botones).
        """
        self._on_select = on_select
        self._close_on_select = close_on_select
        self._title.setText(title.upper())
        # Limpiar celdas anteriores.
        for cell in self._cells:
            cell.setParent(None)
            cell.deleteLater()
        self._cells = []
        self._cell_ids = []
        col_width = (METRICS["flyoutW"] - 22 - 6 * max(columns - 1, 0)) // columns
        for index, item in enumerate(items):
            cell = FlyoutCell(item.icon, item.label, col_width, self)
            cell.setEnabled(item.enabled)
            cell.setProperty("cls", "flyItem")
            if item.tooltip:
                cell.setToolTip(item.tooltip)
            if item.item_id == active_id:
                cell.set_active(True)
            cell.clicked.connect(
                lambda checked=False, iid=item.item_id, trig=item.trigger: self._on_cell(iid, trig)
            )
            self._grid.addWidget(cell, index // columns, index % columns)
            self._cells.append(cell)
            self._cell_ids.append(item.item_id)
        # Pie.
        for btn in self._foot_buttons:
            btn.setParent(None)
            btn.deleteLater()
        self._foot_buttons = []
        if footer is not None:
            self._foot_sep.setVisible(True)
            if footer.text:
                self._foot_txt.setText(footer.text)
                self._foot_txt.setVisible(True)
            else:
                self._foot_txt.setVisible(False)
            for spec in footer.buttons[:3]:
                btn = QToolButton(self)
                btn.setProperty("cls", "flyFoot")
                btn.setText(spec.text)
                btn.setAutoRaise(True)
                btn.setCursor(Qt.CursorShape.PointingHandCursor)
                btn.clicked.connect(
                    lambda checked=False, cb=spec.callback, close=spec.close_after:
                    self._on_footer(cb, close)
                )
                self._foot_row.addWidget(btn)
                self._foot_buttons.append(btn)
        else:
            self._foot_sep.setVisible(False)
            self._foot_txt.setVisible(False)
        self.adjustSize()

    def show_near(self, anchor: QWidget, parent: QWidget) -> None:
        """Muestra el flyout junto a ``anchor`` (a la derecha, clampado)."""
        self.setParent(parent)
        self.show()
        self.update()
        anchor_pos = anchor.mapTo(parent, anchor.rect().topRight())
        x = anchor_pos.x() + 8
        y = anchor_pos.y() - self.height() // 2
        # Clamp al rect de la ventana.
        margin = 8
        max_x = parent.width() - self.width() - margin
        min_x = margin
        x = max(min_x, min(x, max_x))
        max_y = parent.height() - self.height() - margin
        min_y = margin
        y = max(min_y, min(y, max_y))
        self.move(x, y)
        self.raise_()
        # Filtro de cierre (Esc / clic fuera) sobre la ventana.
        self._close_filter_target = parent
        parent.installEventFilter(self)

    def close_with(self, item_id: str | None = None) -> None:
        """Oculta el flyout (emite ``closed``) y opcionalmente notifica la selección."""
        was_visible = self.isVisible()
        self.hide()
        self.remove_close_filter()
        if was_visible:
            self.closed.emit()

    # ------------------------------------------------------------------
    # Internos
    # ------------------------------------------------------------------
    def _on_cell(self, item_id: str, trigger: Callable[[], None] | None = None) -> None:
        if trigger is not None:
            trigger()
        if self._on_select is not None:
            self._on_select(item_id)
        if self._close_on_select:
            self.close_with(item_id)

    def _on_footer(self, callback: Callable[[], None], close: bool) -> None:
        callback()
        if close:
            self.close_with(None)

    def remove_close_filter(self) -> None:
        if self._close_filter_target is not None:
            self._close_filter_target.removeEventFilter(self)
            self._close_filter_target = None

    def keyPressEvent(self, e) -> None:  # noqa: N802
        if e.key() == Qt.Key.Key_Escape and self.isVisible():
            self.close_with(None)
            e.accept()
            return
        super().keyPressEvent(e)

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:
        """Cierra con ``Esc`` o clic fuera (solo mientras es visible)."""
        if not self.isVisible() or obj is not self._close_filter_target:
            return False
        if event.type() == QEvent.Type.KeyPress and event.key() == Qt.Key.Key_Escape:
            self.close_with(None)
            return True
        if event.type() == QEvent.Type.MouseButtonPress:
            me = event  # QMouseEvent
            if me.button() == Qt.MouseButton.LeftButton:
                pos = self.mapFromGlobal(me.globalPos().toPoint())
                if not self.rect().contains(pos):
                    self.close_with(None)
                    return True
        return False
