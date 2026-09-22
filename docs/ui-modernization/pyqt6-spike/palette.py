"""Paleta de comandos (Ctrl+K) para el spike.

Overlay hijo de la ventana principal con fondo translúcido y tarjeta central.
QSS no permite un overlay "modal sobre la propia ventana" sin un QDialog, pero
un QWidget hijo que cubre la ventana y captura teclado (via eventFilter) es
más fiel al mockup (backdrop sobre la app, no un diálogo del SO).

- filtro por substring sobre (sección + título)
- ↑/↓ navega, Enter ejecuta, Esc cierra
- clic en una fila también ejecuta; hover resalta
"""
from __future__ import annotations

from typing import Callable, Sequence

from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QLineEdit, QScrollArea, QVBoxLayout, QWidget,
)

from icons import IconProvider
from theme import ThemeGetter
from widgets import Kbd, shadow


class PaletteAction:
    def __init__(self, section: str, title: str, icon: str, kbd: str = "",
                 fn: Callable[[], None] | None = None):
        self.section = section
        self.title = title
        self.icon = icon
        self.kbd = kbd
        self.fn = fn

    def matches(self, q: str) -> bool:
        return q in f"{self.section} {self.title}".lower()


class CommandPalette(QWidget):
    def __init__(self, icons: IconProvider, theme_getter: ThemeGetter,
                 actions: Sequence[PaletteAction], parent: QWidget):
        super().__init__(parent)
        self.icons, self._tg = icons, theme_getter
        self.actions = list(actions)
        self._filtered: list[PaletteAction] = list(self.actions)
        self._idx = 0
        self._rows: list[tuple[QFrame, int]] = []

        self.setObjectName("paletteOverlay")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._visible = False
        self.hide()

        # --- tarjeta -----------------------------------------------------
        self.card = QFrame(self)
        self.card.setObjectName("paletteCard")
        self.card.setFixedWidth(560)
        shadow(self.card, theme_getter, blur=48, dy=18)
        card_lay = QVBoxLayout(self.card)
        card_lay.setContentsMargins(0, 0, 0, 0)
        card_lay.setSpacing(0)

        # fila de entrada
        self.input_row = QFrame(self.card)
        self.input_row.setObjectName("palInputRow")
        ir = QHBoxLayout(self.input_row)
        ir.setContentsMargins(15, 13, 15, 13)
        ir.setSpacing(10)
        ic = QLabel(); ic.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        ic.setObjectName("palSearchIc")
        self.input = QLineEdit(self.input_row)
        self.input.setObjectName("palInput")
        self.input.setFrame(False)
        self.input.setPlaceholderText("Buscar o ejecutar un comando…")
        self.input.setTextMargins(0, 0, 0, 0)
        ir.addWidget(ic); ir.addWidget(self.input, 1); ir.addWidget(Kbd("Esc"))
        card_lay.addWidget(self.input_row)

        # lista
        self.scroll = QScrollArea(self.card)
        self.scroll.setObjectName("palScroll")
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll.setFixedHeight(330)
        self._list_host = QWidget()
        self._list_host.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self._list_lay = QVBoxLayout(self._list_host)
        self._list_lay.setContentsMargins(6, 6, 6, 6)
        self._list_lay.setSpacing(1)
        self._list_lay.addStretch(1)
        self.scroll.setWidget(self._list_host)
        card_lay.addWidget(self.scroll)

        self.empty = QLabel("Sin resultados")
        self.empty.setObjectName("palEmpty")
        self.empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty.setContentsMargins(20, 24, 20, 24)
        card_lay.addWidget(self.empty)

        # --- señales -------------------------------------------------------
        self.input.textChanged.connect(self._on_query)
        self.input.returnPressed.connect(self._execute_selected)
        self._input_ic = ic
        self.refresh_theme()

    # ------------------------------------------------------------------
    # visibilidad
    def show_overlay(self) -> None:
        self.setGeometry(self.parent().rect())
        self._visible = True
        self.input.clear()
        self._idx = 0
        self._filtered = list(self.actions)
        self._rebuild()
        self.show()
        self.raise_()
        self.input.setFocus()

    def close_overlay(self) -> None:
        self._visible = False
        self.hide()

    def is_visible(self) -> bool:
        return self._visible

    # ------------------------------------------------------------------
    # filtro + render
    def _on_query(self, q: str) -> None:
        qq = q.strip().lower()
        self._filtered = [a for a in self.actions if a.matches(qq)] if qq else list(self.actions)
        self._idx = 0
        self._rebuild()

    def _clear_list(self) -> None:
        for it in self._rows:
            it[0].deleteLater()
        self._rows = []
        # se limpian los QLabel de sección (no están en _rows)
        for it in self._list_host.findChildren(QLabel):
            it.deleteLater()

    def _rebuild(self) -> None:
        self._clear_list()
        self.empty.setText(f'Sin resultados para “{self.input.text()}”')
        self.empty.setVisible(not self._filtered)
        self.card.setVisible(True)
        t = self._tg()
        last_section: str | None = None
        for i, a in enumerate(self._filtered):
            if a.section != last_section:
                sec = QLabel(a.section.upper(), self._list_host)
                sec.setObjectName("palSec")
                sec.setContentsMargins(10, 8, 10, 3)
                self._list_lay.insertWidget(self._list_lay.count() - 1, sec)
                last_section = a.section
            row = QFrame(self._list_host)
            row.setProperty("cls", "palItem")
            row.setProperty("selected", i == self._idx)
            row.setCursor(Qt.CursorShape.PointingHandCursor)
            rl = QHBoxLayout(row)
            rl.setContentsMargins(10, 9, 10, 9)
            rl.setSpacing(11)
            icl = QLabel()
            icl.setPixmap(self.icons.pixmap(a.icon, t["text2"], 17))
            ttl = QLabel(a.title)
            ttl.setObjectName("palTitle")
            rl.addWidget(icl); rl.addWidget(ttl, 1)
            if a.kbd:
                kb = QLabel(a.kbd); kb.setObjectName("palKbd")
                rl.addWidget(kb)
            row.mousePressEvent = lambda e, i=i: self._execute_index(i)
            self._list_lay.insertWidget(self._list_lay.count() - 1, row)
            self._rows.append((row, i))
        self._scroll_selected_into_view()
        self._adjust_height()

    def _adjust_height(self) -> None:
        self.scroll.setVisible(bool(self._filtered))

    def _scroll_selected_into_view(self) -> None:
        if self._rows:
            row, _ = self._rows[self._idx]
            self.scroll.verticalScrollBar().setValue(
                max(0, row.y() - self.scroll.viewport().height() // 2))

    # ------------------------------------------------------------------
    # navegación
    def _move(self, d: int) -> None:
        if not self._filtered:
            return
        self._idx = (self._idx + d) % len(self._filtered)
        for row, i in self._rows:
            on = i == self._idx
            row.setProperty("selected", on)
            st = row.style(); st.unpolish(row); st.polish(row)
        self._scroll_selected_into_view()

    def _execute_index(self, i: int) -> None:
        if 0 <= i < len(self._filtered):
            self._idx = i
            self._execute_selected()

    def _execute_selected(self) -> None:
        if 0 <= self._idx < len(self._filtered):
            a = self._filtered[self._idx]
            self.close_overlay()
            if a.fn:
                a.fn()

    # ------------------------------------------------------------------
    # teclado
    def keyPressEvent(self, e) -> None:  # noqa: N802
        if not self._visible:
            super().keyPressEvent(e)
            return
        k = e.key()
        if k == Qt.Key.Key_Escape:
            e.accept(); self.close_overlay()
        elif k == Qt.Key.Key_Down:
            e.accept(); self._move(1)
        elif k == Qt.Key.Key_Up:
            e.accept(); self._move(-1)
        elif k in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            e.accept(); self._execute_selected()
        else:
            super().keyPressEvent(e)

    def refresh_theme(self) -> None:
        self._input_ic.setPixmap(self.icons.pixmap("search", self._tg()["text3"], 16))
        self._rebuild()
