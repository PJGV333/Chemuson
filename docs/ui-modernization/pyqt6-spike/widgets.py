"""Widgets personalizados pequeños para el spike.

Todo lo que QSS no cubre (indicadores, sublíneas, glifos de atajo, sombras de
flyout/paleta) se resuelve aquí con widgets mínimos y documentados.

Patrones usados:
- `cls`/`sev`/`active` como propiedades dinámicas → selectores QSS por atributo.
- `theme_getter()` que devuelve el Theme *actual* (inmutable por tema) para que
  los widgets pintados a mano se actualicen al cambiar de tema.
"""
from __future__ import annotations

from typing import Callable, Sequence

from PyQt6.QtCore import QSize, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter, QPen
from PyQt6.QtWidgets import (
    QFrame, QGraphicsDropShadowEffect, QGridLayout, QHBoxLayout, QLabel,
    QProgressBar, QScrollArea, QToolButton, QVBoxLayout, QWidget,
)

from icons import IconProvider
from theme import METRICS, Theme, ThemeGetter


# ---------------------------------------------------------------------------
# utilidades
# ---------------------------------------------------------------------------
def shadow(frame: QFrame, theme_getter: ThemeGetter, blur: int = 40, dy: int = 14) -> None:
    """box-shadow aproximado (QSS no lo soporta). Recreable en refresh_theme."""
    eff = QGraphicsDropShadowEffect(frame)
    eff.setBlurRadius(blur)
    eff.setOffset(0, dy)
    frame.setGraphicsEffect(eff)
    frame.shadow_effic_ = eff  # type: ignore[attr-defined]
    _recolor_shadow(frame, theme_getter)


def _recolor_shadow(frame: QFrame, theme_getter: ThemeGetter) -> None:
    eff = getattr(frame, "shadow_effic_", None)
    if eff is not None:
        eff.setColor(theme_getter()["shadow2"])


class Kbd(QLabel):
    """Tecla estilo <kbd> del mockup."""

    def __init__(self, text: str, parent=None):
        super().__init__(text, parent)
        self.setObjectName("kbdK")
        f = self.font(); f.setPixelSize(10); f.setBold(True); self.setFont(f)


class SearchPill(QFrame):
    """Píldora de búsqueda de la app bar (abre la paleta Ctrl+K)."""

    def __init__(self, icons: IconProvider, theme_getter: ThemeGetter, on_click, parent=None):
        super().__init__(parent)
        self.setObjectName("searchPill")
        self.icons, self._tg = icons, theme_getter
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self._on_click = on_click
        lay = QHBoxLayout(self)
        lay.setContentsMargins(10, 6, 10, 6)
        lay.setSpacing(8)
        self.ic = QLabel()
        self.txt = QLabel("Buscar o ejecutar…"); self.txt.setObjectName("searchPillTxt")
        self.kbd = Kbd("Ctrl K")
        lay.addWidget(self.ic); lay.addWidget(self.txt)
        lay.addStretch(1); lay.addWidget(self.kbd)
        self.refresh_theme()

    def mousePressEvent(self, e) -> None:  # noqa: N802
        if e.button() == Qt.MouseButton.LeftButton:
            self._on_click()

    def refresh_theme(self) -> None:
        self.ic.setPixmap(self.icons.pixmap("search", self._tg()["text3"], 14))


# ---------------------------------------------------------------------------
# rail
# ---------------------------------------------------------------------------
class RailButton(QToolButton):
    """Botón del rail (42 px). Icono SVG re-tinted por estado + atajo pintado.

    QSS controla fondo/borde/radio; el color del icono y la letra del atajo se
    pintan a mano porque los SVG están teñidos a color fijo.
    """

    def __init__(self, icons: IconProvider, theme_getter: ThemeGetter, icon_key: str,
                 key_hint: str = "", tooltip: str = "", parent=None):
        super().__init__(parent)
        self.icons, self._tg = icons, theme_getter
        self._icon_key = icon_key
        self._key = key_hint
        self.setProperty("cls", "rail")
        if tooltip:
            self.setToolTip(tooltip)
        self.setFixedSize(METRICS["railBtn"], METRICS["railBtn"])
        self.setIconSize(QSize(METRICS["railIcon"], METRICS["railIcon"]))
        self._refresh_icon()

    def _color(self) -> str:
        if self.property("active"):
            return self._tg()["iconActive"]
        if self.underMouse():
            return self._tg()["iconHover"]
        return self._tg()["icon"]

    def _refresh_icon(self) -> None:
        self.setIcon(self.icons.icon(self._icon_key, self._color(), METRICS["railIcon"]))

    def set_active(self, on: bool) -> None:
        self.setProperty("active", on)
        st = self.style(); st.unpolish(self); st.polish(self)
        self._refresh_icon(); self.update()

    def enterEvent(self, e) -> None:  # noqa: N802
        self._refresh_icon(); super().enterEvent(e)

    def leaveEvent(self, e) -> None:  # noqa: N802
        self._refresh_icon(); super().leaveEvent(e)

    def paintEvent(self, e) -> None:  # noqa: N802
        super().paintEvent(e)
        if not self._key:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        f = self.font(); f.setPixelSize(8); f.setBold(True)
        p.setFont(f)
        color = self._tg()["accentStrong"] if self.property("active") else self._tg()["text3"]
        p.setPen(QPen(QColor(color)))
        p.drawText(QRectF(0, self.height() - 13, self.width() - 4, 11),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, self._key)
        p.end()


# ---------------------------------------------------------------------------
# flyout
# ---------------------------------------------------------------------------
FlyoutItem = tuple[str, str, str]  # (id, label, visual: clave SVG o "glyph:<texto>")


def _wrap_text(text: str, fm, width: int) -> str:
    """Word-wrap greedy por espacios + ruptura de palabra larga.

    El `wordWrap` de QLabel no rompe palabras más anchas que el ancho
    disponible (se desbordan y se recortan); esto replica lo que haría el
    navegador en `.flyout-item .lbl`.
    """
    lines: list[str] = []
    cur: list[str] = []
    cur_w = 0

    def flush() -> None:
        nonlocal cur, cur_w
        if cur:
            lines.append(" ".join(cur))
        cur, cur_w = [], 0

    for word in text.split(" "):
        ww = fm.horizontalAdvance(word)
        if ww > width:
            flush()  # palabra larga: romperla por caracteres
            piece = ""
            for ch in word:
                if fm.horizontalAdvance(piece + ch) <= width:
                    piece += ch
                else:
                    lines.append(piece)
                    piece = ch
            if piece:
                lines.append(piece)
            continue
        trial_w = cur_w + (fm.horizontalAdvance(" ") if cur else 0) + ww
        if trial_w <= width:
            cur.append(word)
            cur_w = trial_w
        else:
            flush()
            cur = [word]
            cur_w = ww
    flush()
    return "\n".join(lines)


class FlyoutCell(QFrame):
    """Celda de flyout: icono 22 px + etiqueta con word-wrap (1-3 líneas).

    Equivalente HTML: `.flyout-item` (flex column, `svg 22px`, `.lbl` 10.5 px
    con `text-align: center` y *sin* `white-space: nowrap` → el texto se
    envuelve). QToolButton no envuelve el texto (elide), por esto es un QFrame
    custom con QLabel word-wrap.
    """
    clicked = pyqtSignal()

    def __init__(self, icons: IconProvider, theme_getter: ThemeGetter,
                 label: str, visual: str, col_w: int, parent=None):
        super().__init__(parent)
        self.icons, self._tg = icons, theme_getter
        self._visual = visual
        self.setProperty("cls", "flyItem")
        self.setToolTip(label)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedWidth(col_w)
        # altura mínima de una línea: sin ella, cuando QGridLayout recalcula
        # mal tras un repopulate, las filas se colapsan a 0 px
        self.setMinimumHeight(56)
        v = QVBoxLayout(self)
        v.setContentsMargins(2, 8, 2, 7)  # HTML: padding 8px 2px 7px
        v.setSpacing(5)                   # HTML: gap 5px
        self.ic = QLabel()
        self.ic.setFixedSize(22, 22)      # HTML: svg 22x22
        self.ic.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl = QLabel(label)
        self.lbl.setObjectName("flyLbl")
        self.lbl.setWordWrap(True)
        self.lbl.setAlignment(
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
        v.addWidget(self.ic, 0, Qt.AlignmentFlag.AlignHCenter)
        v.addWidget(self.lbl)
        # envoltura manual (10 px, como el QSS #flyLbl): el QLabel solo
        # rompe por espacios; las palabras largas se recortarían
        from PyQt6.QtGui import QFontMetrics
        f10 = self.lbl.font(); f10.setPixelSize(10)
        self.lbl.setText(_wrap_text(label, QFontMetrics(f10), col_w - 6))

    def click(self) -> None:
        """API de QToolButton para el smoke test."""
        self.clicked.emit()

    def mousePressEvent(self, e) -> None:  # noqa: N802
        if e.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(e)

    def set_active(self, on: bool) -> None:
        self.setProperty("active", on)
        st = self.style(); st.unpolish(self); st.polish(self)
        t = self._tg()
        color = t["accentStrong"] if on else t["text2"]
        self.ic.setPixmap(self.icons.pixmap(self._visual, color, 22))
        self.update()


class Flyout(QFrame):
    """Flyout de paleta: cabecera + cuadrícula de opciones + pie opcional.

    Ancho fijo 244 px (HTML: `.flyout { width: 244px; padding: 11px }`);
    la cuadrícula divide ese ancho en N columnas como el `grid repeat(N, 1fr)`.
    Se muestra con `show_near(...)`; se cierra con Esc, al elegir una opción o
    al hacer clic fuera (event filter global que instala el app).
    """

    def __init__(self, icons: IconProvider, theme_getter: ThemeGetter, parent=None):
        super().__init__(parent)
        self.icons, self._tg = icons, theme_getter
        self.setObjectName("flyout")
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setFixedWidth(METRICS["flyoutW"])
        self._items: list[QToolButton] = []
        self._ids: list[str] = []
        self._visuals: list[str] = []
        self._active_id: str | None = None
        self._group: str = ""
        self.on_select: Callable[[str, str], None] | None = None
        self.hide()

        lay = QVBoxLayout(self)
        lay.setContentsMargins(11, 11, 11, 11)
        lay.setSpacing(9)

        head = QHBoxLayout(); head.setSpacing(8)
        self.title = QLabel(); self.title.setObjectName("flyoutTitle")
        self._kbd = Kbd("Esc")
        head.addWidget(self.title); head.addStretch(1); head.addWidget(self._kbd)
        lay.addLayout(head)

        self.grid = QGridLayout(); self.grid.setSpacing(6)
        lay.addLayout(self.grid)

        self.foot_sep = QFrame(); self.foot_sep.setObjectName("flyoutSep")
        self.foot_sep.setFixedHeight(1)
        self.foot_txt = QLabel(); self.foot_txt.setObjectName("flyoutFootTxt")
        self.foot_btn = QToolButton()
        self.foot_btn.setProperty("cls", "flyFoot")
        self.foot_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._foot_action: str | None = None
        self.foot_btn.clicked.connect(self._on_foot)
        foot = QHBoxLayout(); foot.setSpacing(8)
        foot.addWidget(self.foot_txt); foot.addStretch(1); foot.addWidget(self.foot_btn)
        self.foot_box = QWidget(); self.foot_box.setLayout(foot)
        lay.addWidget(self.foot_sep); lay.addWidget(self.foot_box)

        shadow(self, theme_getter)

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    def _replace_grid(self) -> QGridLayout:
        """Recrea la cuadrícula en cada populate.

        Reutilizar un QGridLayout tras vaciarlo (takeAt × N) deja
        `totalSize()` obsoleto (0) en Qt 6.11; con el flyout *visible*,
        `adjustSize()` usa ese sizeHint y colapsa el flyout (filas a 0 px).
        Un grid nuevo parte de estado limpio y su sizeHint se calcula bien
        (comportamiento igual al primer populate, que siempre funciona).
        """
        outer = self.layout()
        idx = 0
        for i in range(outer.count()):
            if outer.itemAt(i).layout() is self.grid:
                idx = i
                break
        outer.takeAt(idx)  # quita el QLayoutItem del layout padre
        old = self.grid
        self.grid = QGridLayout()
        self.grid.setSpacing(6)
        outer.insertLayout(idx, self.grid)
        old.deleteLater()
        return self.grid

    def populate(self, title: str, cols: int, items: Sequence[FlyoutItem],
                 foot: tuple[str, str] | None = None,
                 active_id: str | None = None) -> None:
        self.title.setText(title.upper())
        self._active_id = active_id
        for it in self._items:
            it.hide()          # fuera de render de inmediato (evita solapamiento
            it.deleteLater()   # mientras el deferred delete se procesa)
        self._items = []
        self._replace_grid()
        # HTML: width 244 + padding 11 → contenido 222; grid repeat(N, 1fr),
        # gap 6 → ancho de columna = (222 - (N-1)*6) / N
        content_w = METRICS["flyoutW"] - 22
        col_w = (content_w - (cols - 1) * self.grid.spacing()) // cols
        for i, (iid, label, visual) in enumerate(items):
            cell = FlyoutCell(self.icons, self._tg, label, visual, col_w, self)
            cell.clicked.connect(lambda _=False, id=iid, lb=label: self._pick(id, lb))
            self.grid.addWidget(cell, i // cols, i % cols)
            self._items.append(cell)
            self._ids.append(iid)
            self._visuals.append(visual)
        self._foot_action = None
        if foot:
            self.foot_box.show(); self.foot_sep.show()
            self.foot_btn.setText(foot[0])
            self.foot_txt.setText(foot[1] if len(foot) > 1 else "")
        else:
            self.foot_box.hide(); self.foot_sep.hide()
        self._set_content_height(len(items), cols, foot is not None)
        self._restyle_items()

    def _set_content_height(self, n_items: int, cols: int, has_foot: bool) -> None:
        """Altura explícita, calculada de los sizeHints de las celdas.

        Con el flyout *visible y con parent*, el sizeHint del layout
        (QVBoxLayout/QGridLayout) queda inválido (QSize()) hasta que el
        event loop procesa un layout pass; `adjustSize()` usaba ese valor y
        colapsaba el flyout a 41 px. Los sizeHints de las celdas (QFrame +
        QVBoxLayout simples) son siempre válidos, así que la altura se
        calcula directamente — sizing por contenido, igual que el HTML.
        """
        rows = (n_items + cols - 1) // cols
        grid_h = 0
        for r in range(rows):
            row = self._items[r * cols:(r + 1) * cols]
            grid_h += max(c.sizeHint().height() for c in row) + (6 if r else 0)
        head_h = max(self.title.sizeHint().height(), self._kbd.sizeHint().height())
        h = 11 + head_h + 9 + grid_h  # margins 11 + spacing 9 (QVBoxLayout)
        if has_foot:
            foot_h = max(self.foot_txt.sizeHint().height(),
                         self.foot_btn.sizeHint().height())
            h += 9 + 1 + foot_h  # spacing + separador 1 px + contenido
        h += 11
        self.setFixedHeight(h)

    def _pick(self, iid: str, label: str) -> None:
        self._active_id = iid
        if self.on_select:
            self.on_select(iid, label)

    def _on_foot(self) -> None:
        if self.on_select and self._foot_action:
            self.on_select(self._foot_action, "pie")

    def show_near(self, anchor: QWidget, group: str, title: str, cols: int,
                  items: Sequence[FlyoutItem], foot: tuple[str, str] | None = None,
                  active_id: str | None = None) -> None:
        self._group = group
        self.populate(title, cols, items, foot, active_id)
        self.show()
        self.raise_()
        from PyQt6.QtWidgets import QApplication
        # mockup: .flyout { left: 66px } = ancho del rail (58) + 8; el botón
        # empieza en rail_x+8, así que x = a.x() queda exactamente en 58+8
        a = anchor.mapToGlobal(anchor.rect().bottomLeft())
        x = a.x()
        screen = QApplication.primaryScreen().availableGeometry()
        y = max(8, min(a.y() - self.height() // 2, screen.bottom() - self.height() - 8))
        self.move(x, y)

    def _restyle_items(self) -> None:
        for c, iid in zip(self._items, self._ids):
            c.set_active(iid == self._active_id)

    def set_active(self, iid: str | None) -> None:
        self._active_id = iid
        self._restyle_items()

    def keyPressEvent(self, e) -> None:  # noqa: N802
        if e.key() == Qt.Key.Key_Escape:
            self.hide()
        super().keyPressEvent(e)

    def refresh_theme(self) -> None:
        self._restyle_items()
        _recolor_shadow(self, self._tg)
        for child in self.findChildren(QLabel):
            child.style().unpolish(child); child.style().polish(child)


# ---------------------------------------------------------------------------
# panel derecho
# ---------------------------------------------------------------------------
class _SideTabStrip(QFrame):
    """Strip interior: fondo opaco surface + sublínea de acento del tab activo.

    El fondo se pinta a mano (sin cadena de transparencia): una cadena de
    `WA_TranslucentBackground` + viewport transparente no se compone en todas
    las plataformas (offscreen la pinta con #000000/#ffffff sin definir).
    """

    def __init__(self, row: "SideTabRow"):
        super().__init__()
        self._row = row

    def paintEvent(self, e) -> None:  # noqa: N802
        super().paintEvent(e)
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        t = self._row._tg()
        p.fillRect(self.rect(), QColor(t["surface"]))
        if 0 <= self._row._active < len(self._row._buttons):
            b = self._row._buttons[self._row._active]
            r = b.geometry()
            p.setPen(QPen(Qt.PenStyle.NoPen))
            p.setBrush(QColor(t["accent"]))
            p.drawRoundedRect(QRectF(r.left() + 8, self.height() - 4.5,
                                     max(8.0, r.width() - 16), 2.5), 1.25, 1.25)
        p.end()


class SideTabRow(QFrame):
    """Fila de tabs del panel derecho (sublínea de acento pintada a mano).

    Equivalente HTML: .side-tabs { overflow-x: auto; scrollbar-width: none }.
    Si los tabs no caben en 324 px, la fila se desplaza horizontalmente
    (sin barra: con rueda horizontal / shift+rueda), igual que el mockup.
    """

    def __init__(self, theme_getter: ThemeGetter, parent=None):
        super().__init__(parent)
        self._tg = theme_getter
        self.setFixedWidth(METRICS["sideW"])
        self.setFixedHeight(40)
        self._buttons: list[QToolButton] = []
        self._active = 0
        self.on_change: Callable[[int], None] | None = None

        self._strip = _SideTabStrip(self)
        self._lay = QHBoxLayout(self._strip)
        self._lay.setContentsMargins(8, 0, 10, 0)
        self._lay.setSpacing(0)
        self._lay.addStretch(1)

        self._sa = QScrollArea(self)
        self._sa.setObjectName("sideTabsScroll")
        self._sa.setFrameShape(QFrame.Shape.NoFrame)
        self._sa.setWidgetResizable(True)
        self._sa.setWidget(self._strip)
        self._sa.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._sa.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(self._sa)

    def add_tab(self, label: str) -> None:
        b = QToolButton(self._strip)
        b.setProperty("cls", "sideTab")
        b.setText(label)
        b.setCheckable(True)
        b.setCursor(Qt.CursorShape.PointingHandCursor)
        b.clicked.connect(lambda _=False, i=len(self._buttons): self._select(i))
        from PyQt6.QtGui import QFontMetrics
        fm = QFontMetrics(b.font())
        b.setMinimumWidth(fm.horizontalAdvance(label) + 20)
        b.setMaximumWidth(10**6)
        self._buttons.append(b)
        self._lay.insertWidget(self._lay.count() - 1, b)

    def _select(self, i: int) -> None:
        self._active = i
        for j, b in enumerate(self._buttons):
            b.setProperty("active", j == i)
            st = b.style(); st.unpolish(b); st.polish(b)
        self._strip.update()
        # mantener visible el tab seleccionado
        if 0 <= i < len(self._buttons):
            b = self._buttons[i]
            sb = self._sa.horizontalScrollBar()
            target = b.x() - 8
            sb.setValue(max(0, min(target, sb.maximum() - 40)))
        if self.on_change:
            self.on_change(i)

    def set_active(self, i: int) -> None:
        self._select(i)

    def _update_tab_widths(self) -> None:
        """Recalcula el ancho de cada tab con la fuente *actual*.

        Los tabs se crean antes de aplicar la fuente del tema; sin esto el
        minimumWidth queda corto y QToolButton elide la etiqueta («Ins...tor»).
        """
        from PyQt6.QtGui import QFontMetrics
        for b in self._buttons:
            fm = QFontMetrics(b.font())
            b.setMinimumWidth(fm.horizontalAdvance(b.text()) + 20)

    def refresh_theme(self) -> None:
        self._update_tab_widths()
        for b in self._buttons:
            st = b.style(); st.unpolish(b); st.polish(b)
        self.update()
        self._strip.update()

    def paintEvent(self, e) -> None:  # noqa: N802
        super().paintEvent(e)
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        t = self._tg()
        p.setPen(QPen(QColor(t["border"])))
        p.drawLine(0, self.height() - 1, self.width(), self.height() - 1)
        p.end()

    def wheelEvent(self, e) -> None:  # noqa: N802
        # rueda vertical → scroll horizontal de los tabs (comodidad; el HTML
        # usa overflow-x auto con scrollbar oculta)
        sb = self._sa.horizontalScrollBar()
        if sb.maximum() > 0:
            sb.setValue(sb.value() + (12 if e.angleDelta().y() < 0 else -12))
            e.accept()
        else:
            super().wheelEvent(e)


# ---------------------------------------------------------------------------
# contenido de paneles
# ---------------------------------------------------------------------------
def make_pill(text: str, sev: str = "") -> QFrame:
    f = QFrame(); f.setProperty("cls", "pill")
    if sev:
        f.setProperty("sev", sev)
    f.setFrameShape(QFrame.Shape.NoFrame)
    lay = QHBoxLayout(f); lay.setContentsMargins(9, 2, 9, 2); lay.setSpacing(0)
    lbl = QLabel(text)
    lf = lbl.font(); lf.setPixelSize(10); lf.setBold(True); lbl.setFont(lf)
    lay.addWidget(lbl)
    return f


def make_kv(k: str, v: str, small: str = "") -> QFrame:
    f = QFrame(); f.setProperty("cls", "kv"); f.setFrameShape(QFrame.Shape.NoFrame)
    lay = QVBoxLayout(f); lay.setContentsMargins(10, 8, 10, 8); lay.setSpacing(2)
    a = QLabel(k); a.setObjectName("kvK")
    b = QLabel(v); b.setObjectName("kvV")
    lay.addWidget(a); lay.addWidget(b)
    if small:
        c = QLabel(small); c.setObjectName("kvSmall")
        lay.addWidget(c)
    return f


def make_row(k: str, right: QWidget | str, mono: bool = False) -> QFrame:
    f = QFrame(); f.setProperty("cls", "row"); f.setFrameShape(QFrame.Shape.NoFrame)
    lay = QHBoxLayout(f); lay.setContentsMargins(2, 7, 2, 7); lay.setSpacing(10)
    a = QLabel(k); a.setObjectName("rowK")
    if isinstance(right, str):
        b = QLabel(right); b.setObjectName("rowV")
        if mono:
            b.setProperty("mono", True)
    else:
        b = right
    lay.addWidget(a); lay.addStretch(1); lay.addWidget(b)
    return f


class Issue(QFrame):
    """Fila de resultado de validación (clicable)."""

    clicked = pyqtSignal()

    def mousePressEvent(self, e) -> None:  # noqa: N802
        if e.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(e)


def make_issue(sev: str, ttl: str, sub: str, chip_text: str,
               theme_getter: ThemeGetter, selected: bool = False) -> Issue:
    f = Issue(); f.setProperty("cls", "issue"); f.setFrameShape(QFrame.Shape.NoFrame)
    f.setCursor(Qt.CursorShape.PointingHandCursor)
    f.setProperty("selected", selected)
    f._sev = sev
    lay = QHBoxLayout(f); lay.setContentsMargins(10, 9, 10, 9); lay.setSpacing(9)
    dot = QFrame(); dot.setFixedSize(8, 8)
    lay.addWidget(dot, 0, Qt.AlignmentFlag.AlignTop)
    body = QWidget(); body.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
    bl = QVBoxLayout(body); bl.setContentsMargins(0, 0, 0, 0); bl.setSpacing(1)
    t = QLabel(ttl); t.setObjectName("issueTtl")
    s = QLabel(sub); s.setObjectName("issueSub")
    bl.addWidget(t); bl.addWidget(s)
    lay.addWidget(body, 1)
    lay.addWidget(make_pill(chip_text, sev), 0, Qt.AlignmentFlag.AlignTop)
    f.dot = dot
    f.refresh_theme = lambda: _recolor_issue(f, theme_getter)  # type: ignore[method-assign]
    _recolor_issue(f, theme_getter)
    return f


def _recolor_issue(f: QFrame, theme_getter: ThemeGetter) -> None:
    t = theme_getter()
    color = {"warn": t["warn"], "err": t["danger"], "ok": t["ok"]}.get(f._sev, t["text3"])
    f.dot.setStyleSheet(f"QFrame {{ background: {color}; border-radius: 4px; }}")


def make_bar(frac: float) -> QProgressBar:
    p = QProgressBar()
    p.setRange(0, 100)
    p.setValue(int(frac * 100))
    p.setTextVisible(False)
    p.setFixedHeight(5)
    p.setMinimumWidth(90)
    return p


def make_tpl_card(icons: IconProvider, theme_getter: ThemeGetter, icon_key: str,
                  name: str, meta: str, on_click=None) -> QFrame:
    f = QFrame(); f.setProperty("cls", "tpl"); f.setFrameShape(QFrame.Shape.NoFrame)
    f.setCursor(Qt.CursorShape.PointingHandCursor)
    lay = QVBoxLayout(f); lay.setContentsMargins(10, 10, 10, 10); lay.setSpacing(7)
    pic = QFrame(); pic.setFixedHeight(46); pic.setAttribute(
        Qt.WidgetAttribute.WA_TranslucentBackground, True)
    pl = QHBoxLayout(pic); pl.setContentsMargins(0, 0, 0, 0); pl.addStretch(1)
    ic = QLabel(); ic.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
    pl.addWidget(ic, 0, Qt.AlignmentFlag.AlignCenter)
    nm = QLabel(name); nm.setObjectName("tplNm")
    mt = QLabel(meta); mt.setObjectName("tplMeta")
    lay.addWidget(pic); lay.addWidget(nm); lay.addWidget(mt)

    def _refresh() -> None:
        ic.setPixmap(icons.pixmap(icon_key, theme_getter()["text2"], 34))

    f.refresh_theme = _refresh  # type: ignore[method-assign]
    _refresh()
    if on_click:
        f.mousePressEvent = lambda e, o=on_click: o()  # type: ignore[method-assign]
    return f


# ---------------------------------------------------------------------------
# chips con punto (Rejilla / Números)
# ---------------------------------------------------------------------------
class ChipToggle(QToolButton):
    """Chip con punto de estado a la derecha (pintado a mano)."""

    def __init__(self, label: str, icons: IconProvider, theme_getter: ThemeGetter,
                 icon_key: str, checked: bool = False, parent=None):
        super().__init__(parent)
        self.icons, self._tg = icons, theme_getter
        self._icon_key = icon_key
        self.setProperty("cls", "chip")
        self.setCheckable(True)
        self.setChecked(checked)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.setIconSize(QSize(13, 13))
        self.setText(f"  {label}  ")
        f = self.font(); f.setPixelSize(11); f.setBold(True); self.setFont(f)
        self.refresh_theme()

    def refresh_theme(self) -> None:
        color = self._tg()["accentStrong"] if self.isChecked() else self._tg()["text2"]
        self.setIcon(self.icons.icon(self._icon_key, color, 13))
        self.update()

    def paintEvent(self, e) -> None:  # noqa: N802
        super().paintEvent(e)
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        color = self._tg()["accent"] if self.isChecked() else self._tg()["text3"]
        p.setPen(QPen(Qt.PenStyle.NoPen))
        p.setBrush(QColor(color))
        p.drawEllipse(QRectF(self.width() - 15, self.height() / 2 - 3, 6, 6))
        p.end()
