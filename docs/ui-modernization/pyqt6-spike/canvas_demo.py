"""Lienzo de demostración: QGraphicsView/QGraphicsScene (misma tecnología que M09).

Reproduce la vista central del mockup:
- fondo de trabajo punteado (estático, en coords de viewport)
- hoja blanca con sombra y radios 4
- rejilla fina dentro de la hoja (en coords de escena: hace zoom con ella)
- molécula de demostración (1-metil-4-nitrobenzeno, coords idénticas al SVG del
  mockup), halos de texto en heteroátomos, numeración, anotación
- rectángulo de selección discontinuo + 4 handles
- pill de zoom y chips Rejilla/Números (en app.py, esta clase expone la API)

No hay lógica química: solo aspecto gráfico.
"""
from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import (
    QBrush, QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen,
)
from PyQt6.QtWidgets import QGraphicsDropShadowEffect, QGraphicsItem, QGraphicsPathItem, QGraphicsRectItem, QGraphicsScene, QGraphicsView

_AA = QPainter.RenderHint.Antialiasing

INK = QColor("#1F2937")          # tinta de la molécula (fija en ambos temas, igual que el mockup)
O_COLOR = QColor("#DC2626")      # oxígeno
N_COLOR = QColor("#2563EB")      # nitrógeno
NUM_COLOR = QColor("#94A3B8")    # numeración
ANNOT = QColor("#0E7490")        # anotación

SHEET_W, SHEET_H = 560, 480
ZOOM_STEPS = (0.5, 0.7, 0.85, 1.0, 1.15, 1.35, 1.6, 2.0)


class AtomLabel(QGraphicsItem):
    """Texto con halo blanco (equivalente CSS: paint-order: stroke).

    QGraphicsItem propio que pinta con `drawText` en dos pasadas (halo blanco
    grueso + color). NOTA: `QPainterPath.addText` no dibuja el texto en este
    binding (el glifo no se convierte en subpath y DrawTextMode no está
    expuesto), por eso el item propio.
    """

    def __init__(self, text: str, x: int, y: int, *, fill: QColor, font_px: int = 19,
                 weight: int = 600, halo_px: float = 5.0, anchor_center_x: bool = True):
        super().__init__()
        self._text = text
        self._fill = fill
        self._halo = halo_px
        self._f = QFont()
        self._f.setPixelSize(font_px)
        self._f.setWeight(weight)
        advance = QFontMetrics(self._f).horizontalAdvance(text)
        tx = x - advance / 2 if anchor_center_x else x
        self._pen_pos = QPointF(tx, y)
        ascent = font_px * 1.4
        self._rect = QRectF(tx - advance / 2 - halo_px, y - ascent - halo_px,
                            advance + 2 * halo_px, ascent + font_px + 2 * halo_px)

    def boundingRect(self) -> QRectF:  # noqa: N802
        return self._rect

    def paint(self, painter: QPainter, option, widget=None) -> None:  # noqa: N802
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        painter.setFont(self._f)
        if self._halo > 0:
            halo = QPen(QColor("#FFFFFF"), self._halo, Qt.PenStyle.SolidLine,
                        Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
            painter.setPen(halo)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawText(self._pen_pos, self._text)
        painter.setPen(QPen(self._fill))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawText(self._pen_pos, self._text)


class MoleculeLabel(QGraphicsItem):
    """Texto SIN halo (numeración pequeña)."""

    def __init__(self, text: str, x: int, y: int, *, fill: QColor = NUM_COLOR, font_px: int = 11):
        super().__init__()
        self._text = text
        self._fill = fill
        self._halo = 0.0
        self._f = QFont()
        self._f.setPixelSize(font_px)
        self._f.setWeight(600)
        advance = QFontMetrics(self._f).horizontalAdvance(text)
        self._pen_pos = QPointF(x - advance / 2, y)
        ascent = font_px * 1.4
        self._rect = QRectF(x - advance / 2, y - ascent, advance, ascent + font_px)

    def boundingRect(self) -> QRectF:  # noqa: N802
        return self._rect

    def paint(self, painter: QPainter, option, widget=None) -> None:  # noqa: N802
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        painter.setFont(self._f)
        painter.setPen(QPen(self._fill))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawText(self._pen_pos, self._text)


def _bond_pen(color: QColor = INK, width: float = 2.4) -> QPen:
    return QPen(color, width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap,
                Qt.PenJoinStyle.RoundJoin)


class CanvasDemo:
    """Escena + vista de demostración."""

    def __init__(self, theme, icon_provider):
        self.theme = theme
        self.icons = icon_provider
        self.scene = QGraphicsScene()
        self.scene.setSceneRect(QRectF(-140, -140, SHEET_W + 280, SHEET_H + 280))
        self.grid_on = True
        self.view = CanvasDemoView(self)
        self._build()

    # ------------------------------------------------------------------
    def _build(self) -> None:
        s = self.scene
        t = self.theme

        # --- hoja ------------------------------------------------------
        # (QGraphicsRectItem.setRadius no está expuesto en este binding de
        #  PyQt6; un QGraphicsPathItem con rect redondeado es equivalente)
        sheet_path = QPainterPath()
        sheet_path.addRoundedRect(QRectF(0, 0, SHEET_W, SHEET_H), 4, 4)
        from PyQt6.QtWidgets import QGraphicsPathItem
        self.sheet = QGraphicsPathItem(sheet_path)
        self.sheet.setBrush(QBrush(QColor(t["sheet"])))
        self.sheet.setPen(QPen(Qt.PenStyle.NoPen))
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(34)
        shadow.setOffset(0, 14)
        shadow.setColor(t["sheetShadow"])
        self.sheet.setGraphicsEffect(shadow)
        self._shadow = shadow
        s.addItem(self.sheet)

        # --- molécula: anillo bencénico + sustituyentes ----------------
        p = QPainterPath()
        # hexágono (coords del mockup)
        ring = [(290, 150), (364, 193), (364, 277), (290, 320), (216, 277), (216, 193)]
        p.moveTo(*ring[0])
        for pt in ring[1:]:
            p.lineTo(*pt)
        p.closeSubpath()
        # círculo aromático
        p.addEllipse(QPointF(290, 235), 44, 44)
        self.ring_item = self._path(p, _bond_pen())
        s.addItem(self.ring_item)

        # enlace C–CH3 + etiqueta
        p = QPainterPath(); p.moveTo(364, 193); p.lineTo(404, 168)
        s.addItem(self._path(p, _bond_pen()))
        s.addItem(AtomLabel("CH\u2083", 418, 172, fill=INK, anchor_center_x=False))

        # grupo NO2: C–N, N=O (doble), N–O
        p = QPainterPath(); p.moveTo(216, 277); p.lineTo(174, 300)
        s.addItem(self._path(p, _bond_pen()))
        p = QPainterPath(); p.moveTo(147, 303); p.lineTo(121, 289)
        p.moveTo(150, 309); p.lineTo(124, 295)
        s.addItem(self._path(p, _bond_pen()))
        p = QPainterPath(); p.moveTo(150, 318); p.lineTo(130, 348)
        s.addItem(self._path(p, _bond_pen()))
        s.addItem(AtomLabel("N", 155, 312, fill=N_COLOR))
        s.addItem(AtomLabel("O", 106, 292, fill=O_COLOR))
        s.addItem(AtomLabel("O", 121, 364, fill=O_COLOR))

        # numeración de átomos (toggleable con el chip Números)
        self.num_items: list[MoleculeLabel] = []
        for num, (x, y) in {
            "1": (290, 140), "2": (376, 188), "3": (376, 284),
            "4": (290, 338), "5": (204, 262), "6": (204, 188),
        }.items():
            lbl = MoleculeLabel(num, x, y)
            s.addItem(lbl)
            self.num_items.append(lbl)

        # --- anotación (flecha curva) ----------------------------------
        p = QPainterPath(); p.moveTo(400, 150); p.cubicTo(380, 118, 356, 112, 302, 136)
        self.annot = self._path(p, QPen(ANNOT, 2.0, Qt.PenStyle.SolidLine,
                                        Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
        s.addItem(self.annot)
        head = QPainterPath(); head.moveTo(302, 136); head.lineTo(313, 127); head.lineTo(315, 139)
        head.closeSubpath()
        self.annot_head = self._path(head, QPen(Qt.PenStyle.NoPen))
        self.annot_head.setBrush(QBrush(ANNOT))
        s.addItem(self.annot_head)

        # --- selección ---------------------------------------------------
        sel = QColor(t["accent"])
        self.sel_rect = QGraphicsRectItem(78, 262, 112, 118)
        s.addItem(self.sel_rect)
        self.sel_rect.setZValue(10)
        self.handles = []
        for hx, hy in ((75, 259), (187, 259), (75, 377), (187, 377)):
            h = QGraphicsRectItem(hx, hy, 7, 7)
            h.setBrush(QBrush(QColor("#FFFFFF")))
            h.setPen(QPen(sel, 1.5))
            h.setZValue(11)
            s.addItem(h)
            self.handles.append(h)
        self._apply_selection_style()

        # --- vista -------------------------------------------------------
        self.view.setRenderHints(_AA | QPainter.RenderHint.SmoothPixmapTransform)
        from PyQt6.QtWidgets import QFrame as _QF; self.view.setFrameShape(_QF.Shape.NoFrame)
        self._fit()

    def _path(self, path: QPainterPath, pen: QPen):
        item = QGraphicsPathItem(path)
        item.setPen(pen)
        return item

    # ------------------------------------------------------------------
    def _apply_selection_style(self) -> None:
        sel = QColor(self.theme["accent"])
        pen = QPen(sel, 1.6, Qt.PenStyle.DashLine)
        pen.setDashPattern([5.0, 4.0])
        self.sel_rect.setPen(pen)
        self.sel_rect.setBrush(QBrush(QColor(14, 116, 144, 13)))
        for h in self.handles:
            h.setPen(QPen(sel, 1.5))

    def _fit(self) -> None:
        self._zoom_idx = 3
        self.view.resetTransform()
        self.view.centerOn(QPointF(SHEET_W / 2, SHEET_H / 2))

    # ------------------------------------------------------------------
    # zoom
    def set_zoom_idx(self, idx: int) -> int:
        self._zoom_idx = max(0, min(len(ZOOM_STEPS) - 1, idx))
        self.view.resetTransform()
        self.view.scale(ZOOM_STEPS[self._zoom_idx], ZOOM_STEPS[self._zoom_idx])
        self.view.centerOn(QPointF(SHEET_W / 2, SHEET_H / 2))
        return self._zoom_idx

    def zoom_in(self) -> int:
        return self.set_zoom_idx(getattr(self, "_zoom_idx", 3) + 1)

    def zoom_out(self) -> int:
        return self.set_zoom_idx(getattr(self, "_zoom_idx", 3) - 1)

    def zoom_fit(self) -> int:
        return self.set_zoom_idx(3)

    def zoom_label(self) -> str:
        return f"{round(ZOOM_STEPS[getattr(self, '_zoom_idx', 3)] * 100)} %"

    # ------------------------------------------------------------------
    # tema
    def apply_theme(self, theme) -> None:
        self.theme = theme
        t = theme.tokens
        self.scene.setBackgroundBrush(QBrush(QColor(t["canvasBg"])))
        self._shadow.setColor(t["sheetShadow"])
        self._apply_selection_style()
        self.view.grid_on = self.grid_on
        self.view.viewport().update()

    def set_grid(self, on: bool) -> None:
        self.grid_on = on
        self.view.grid_on = on
        self.view.viewport().update()

    def set_numbers(self, on: bool) -> None:
        for it in self.num_items:
            it.setVisible(on)


class CanvasDemoView(QGraphicsView):
    """Vista que pinta el fondo punteado (viewport) y la rejilla de la hoja (escena)."""

    cursor_moved = pyqtSignal(QPointF)

    def __init__(self, canvas: CanvasDemo):
        super().__init__(canvas.scene)
        self._canvas = canvas
        self.grid_on = True
        self.setObjectName("canvas")
        self.setRenderHint(_AA)

    # ------------------------------------------------------------------
    def drawBackground(self, painter, rect: QRectF) -> None:  # noqa: N802 (API Qt)
        t = self._canvas.theme.tokens
        # 1) fondo + puntos: coordenadas de VIEWPORT (estático, como el CSS del mockup)
        painter.save()
        painter.resetTransform()
        painter.fillRect(painter.viewport(), QColor(t["canvasBg"]))
        # los puntos del fondo de trabajo son independientes de la rejilla de la hoja
        dot = QColor(t["border"])
        painter.setPen(QPen(dot, 1.6))
        vp_w, vp_h = painter.viewport().width(), painter.viewport().height()
        step = 22
        x = 11
        while x < vp_w:
            y = 11
            while y < vp_h:
                painter.drawPoint(int(x), int(y))
                y += step
            x += step
        painter.restore()
        # 2) rejilla de la hoja: coordenadas de ESCENA (hace zoom/pan con la hoja)
        if self._canvas.grid_on:
            painter.save()
            sheet = QRectF(0, 0, SHEET_W, SHEET_H)
            painter.setClipRect(sheet)
            grid = QColor(t["sheetGrid"])
            painter.setPen(QPen(grid, 1.0))
            step = 24
            x0, x1 = int(sheet.left()), int(sheet.right())
            y0, y1 = int(sheet.top()), int(sheet.bottom())
            for x in range(x0, x1 + 1, step):
                painter.drawLine(QPointF(x, y0), QPointF(x, y1))
            for y in range(y0, y1 + 1, step):
                painter.drawLine(QPointF(x0, y), QPointF(x1, y))
            painter.restore()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        super().mouseMoveEvent(event)
        self.cursor_moved.emit(self.mapToScene(event.position().toPoint()))
