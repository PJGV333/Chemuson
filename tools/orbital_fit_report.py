"""Reporte geométrico de ajuste orbital contra una hoja de referencia 2D."""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QApplication

from chemuson.gui.orbitals import (
    ORBITAL_PALETTE_MODEL,
    ORBITAL_SPECS,
    GlyphLayer,
    orbital_renderer,
)


REFERENCE_SHEET = ROOT / "tests" / "data" / "orbitals" / "reference_sheet.png"
OUTPUT_DIR = ROOT / "tests" / "data" / "orbitals" / "fit_report"
REFERENCE_CROPS_DIR = OUTPUT_DIR / "reference_crops"
CURRENT_CROPS_DIR = OUTPUT_DIR / "current_crops"
OVERLAY_DIR = OUTPUT_DIR / "overlay"
REPORT_PATH = OUTPUT_DIR / "report.txt"
GRID_PREVIEW_PATH = OUTPUT_DIR / "grid_preview.png"
OVERLAY_SHEET_PATH = OUTPUT_DIR / "overlay_sheet.png"

CANVAS_SIZE = 96
BG = QColor("#FFFFFF")


def _ellipse(x: float, y: float, w: float, h: float) -> QPainterPath:
    path = QPainterPath()
    path.addEllipse(QRectF(x, y, w, h))
    return path


def _ring(x: float, y: float, w: float, h: float, inner_x: float, inner_y: float) -> QPainterPath:
    path = QPainterPath()
    path.setFillRule(Qt.FillRule.OddEvenFill)
    path.addEllipse(QRectF(x, y, w, h))
    inset_x = w * (1.0 - inner_x) * 0.5
    inset_y = h * (1.0 - inner_y) * 0.5
    path.addEllipse(QRectF(x + inset_x, y + inset_y, w * inner_x, h * inner_y))
    return path


def _cubic(
    start: tuple[float, float],
    segments: tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], ...],
) -> QPainterPath:
    path = QPainterPath(QPointF(*start))
    for c1, c2, end in segments:
        path.cubicTo(QPointF(*c1), QPointF(*c2), QPointF(*end))
    path.closeSubpath()
    return path


def _transform(path: QPainterPath, *, rotate: float = 0.0, tx: float = 0.0, ty: float = 0.0, sx: float = 1.0, sy: float = 1.0) -> QPainterPath:
    from PyQt6.QtGui import QTransform

    transform = QTransform()
    transform.translate(tx, ty)
    if rotate:
        transform.rotate(rotate)
    transform.scale(sx, sy)
    return transform.map(path)


def _draw_layers(image: QImage, layers: tuple[GlyphLayer, ...]) -> QImage:
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    for layer in layers:
        if layer.paint == "outline":
            painter.setPen(QPen(QColor("#111111"), 2.1, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
            painter.setBrush(Qt.BrushStyle.NoBrush)
        elif layer.paint == "shaded":
            if layer.stroke:
                painter.setPen(QPen(QColor("#111111"), 1.7, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
            else:
                painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor("#7A7A7A"))
        else:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor("#111111"))
        painter.drawPath(layer.path)
    painter.end()
    return image


def _reference_layers(kind: str) -> tuple[GlyphLayer, ...]:
    cx = CANVAS_SIZE * 0.5
    cy = CANVAS_SIZE * 0.5

    def move(path: QPainterPath) -> QPainterPath:
        return _transform(path, tx=cx, ty=cy)

    p_top = _cubic(
        (0.0, -41.0),
        (
            ((10.8, -39.5), (15.2, -29.0), (13.6, -19.5)),
            ((11.2, -9.0), (5.4, -2.6), (0.0, 0.0)),
            ((-5.4, -2.6), (-11.2, -9.0), (-13.6, -19.5)),
            ((-15.2, -29.0), (-10.8, -39.5), (0.0, -41.0)),
        ),
    )
    p_bottom = _transform(p_top, sy=-1.0)

    sp_small = _cubic(
        (0.0, -26.0),
        (
            ((3.8, -25.0), (5.3, -20.5), (4.4, -15.8)),
            ((3.4, -11.2), (1.7, -8.2), (0.0, -6.0)),
            ((-1.7, -8.2), (-3.4, -11.2), (-4.4, -15.8)),
            ((-5.3, -20.5), (-3.8, -25.0), (0.0, -26.0)),
        ),
    )
    sp_lobe = _cubic(
        (0.0, -6.0),
        (
            ((11.8, -4.6), (15.0, 8.2), (12.4, 22.0)),
            ((10.4, 31.5), (5.2, 37.2), (0.0, 41.0)),
            ((-5.2, 37.2), (-10.4, 31.5), (-12.4, 22.0)),
            ((-15.0, 8.2), (-11.8, -4.6), (0.0, -6.0)),
        ),
    )

    d_top = _cubic(
        (0.0, -38.0),
        (
            ((10.0, -36.2), (14.8, -25.5), (13.0, -15.6)),
            ((10.8, -7.2), (5.2, -1.8), (0.0, 0.4)),
            ((-5.2, -1.8), (-10.8, -7.2), (-13.0, -15.6)),
            ((-14.8, -25.5), (-10.0, -36.2), (0.0, -38.0)),
        ),
    )
    d_bottom = _transform(d_top, sy=-1.0)
    d_left = _cubic(
        (-36.0, 0.0),
        (
            ((-33.5, -6.2), (-21.0, -10.2), (-8.8, -8.2)),
            ((-2.8, -6.6), (-0.6, -2.8), (0.0, 0.0)),
            ((-0.6, 2.8), (-2.8, 6.6), (-8.8, 8.2)),
            ((-21.0, 10.2), (-33.5, 6.2), (-36.0, 0.0)),
        ),
    )
    d_right = _transform(d_left, sx=-1.0)

    dz2_top = _cubic(
        (0.0, -39.0),
        (
            ((8.2, -37.4), (12.2, -26.0), (10.4, -17.0)),
            ((8.8, -10.0), (4.2, -6.2), (0.0, -5.6)),
            ((-4.2, -6.2), (-8.8, -10.0), (-10.4, -17.0)),
            ((-12.2, -26.0), (-8.2, -37.4), (0.0, -39.0)),
        ),
    )
    dz2_bottom = _transform(dz2_top, sy=-1.0)
    dz2_ring = _ring(-31.0, -4.8, 62.0, 9.6, 0.46, 0.26)

    sigma_large = _cubic(
        (0.0, -38.0),
        (
            ((12.0, -35.0), (16.0, -20.0), (14.0, -8.0)),
            ((12.0, 3.0), (7.0, 15.0), (0.0, 28.0)),
            ((-7.0, 15.0), (-12.0, 3.0), (-14.0, -8.0)),
            ((-16.0, -20.0), (-12.0, -35.0), (0.0, -38.0)),
        ),
    )
    sigma_small = _cubic(
        (0.0, 22.0),
        (
            ((4.5, 21.0), (6.5, 24.0), (6.0, 28.0)),
            ((5.0, 31.0), (2.5, 33.0), (0.0, 34.0)),
            ((-2.5, 33.0), (-5.0, 31.0), (-6.0, 28.0)),
            ((-6.5, 24.0), (-4.5, 21.0), (0.0, 22.0)),
        ),
    )

    pi_top = _transform(dz2_top, ty=-2.0)
    pi_bottom = _transform(pi_top, sy=-1.0)
    pi_ring = _ring(-24.0, -5.5, 48.0, 11.0, 0.44, 0.40)

    torus_outline = (GlyphLayer(move(_ring(-36.0, -10.0, 72.0, 20.0, 0.40, 0.34)), "outline"),)
    torus_shaded = (GlyphLayer(move(_ring(-36.0, -10.0, 72.0, 20.0, 0.40, 0.34)), "shaded"),)
    torus_solid = (GlyphLayer(move(_ring(-36.0, -10.0, 72.0, 20.0, 0.40, 0.34)), "solid"),)

    reference = {
        "p_shaded": (
            GlyphLayer(move(p_top), "shaded", stroke=False),
            GlyphLayer(move(p_bottom), "outline"),
        ),
        "p_solid": (GlyphLayer(move(p_top), "solid"), GlyphLayer(move(p_bottom), "outline")),
        "sp_lobe_outline": (GlyphLayer(move(sp_lobe), "outline"),),
        "sp_lobe_shaded": (
            GlyphLayer(move(sp_small), "shaded", stroke=False),
            GlyphLayer(move(sp_lobe), "shaded", stroke=False),
        ),
        "sp_lobe_solid": (
            GlyphLayer(move(sp_small), "solid"),
            GlyphLayer(move(sp_lobe), "solid"),
        ),
        "d_shaded": (
            GlyphLayer(move(d_left), "outline"),
            GlyphLayer(move(d_top), "shaded", stroke=False),
            GlyphLayer(move(d_right), "outline"),
            GlyphLayer(move(d_bottom), "shaded", stroke=False),
        ),
        "d_solid": (
            GlyphLayer(move(d_left), "outline"),
            GlyphLayer(move(d_top), "solid"),
            GlyphLayer(move(d_right), "outline"),
            GlyphLayer(move(d_bottom), "solid"),
        ),
        "dz2_shaded": (
            GlyphLayer(move(dz2_top), "shaded", stroke=False),
            GlyphLayer(move(dz2_ring), "shaded", stroke=False),
            GlyphLayer(move(dz2_bottom), "shaded", stroke=False),
        ),
        "dz2_solid": (
            GlyphLayer(move(dz2_top), "solid"),
            GlyphLayer(move(dz2_ring), "solid"),
            GlyphLayer(move(dz2_bottom), "solid"),
        ),
        "sigma_bonding_shaded": (
            GlyphLayer(move(sigma_large), "outline"),
            GlyphLayer(move(sigma_small), "shaded"),
        ),
        "sigma_bonding_solid": (
            GlyphLayer(move(sigma_large), "outline"),
            GlyphLayer(move(sigma_small), "solid"),
        ),
        "pi_bonding_shaded": (
            GlyphLayer(move(pi_top), "outline"),
            GlyphLayer(move(pi_ring), "shaded"),
            GlyphLayer(move(pi_bottom), "outline"),
        ),
        "pi_bonding_solid": (
            GlyphLayer(move(pi_top), "outline"),
            GlyphLayer(move(pi_ring), "solid"),
            GlyphLayer(move(pi_bottom), "outline"),
        ),
        "torus_outline": torus_outline,
        "torus_shaded": torus_shaded,
        "torus_solid": torus_solid,
    }
    return reference.get(kind, ())


def _render_reference_crop(kind: str) -> QImage:
    image = QImage(CANVAS_SIZE, CANVAS_SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(BG)
    layers = _reference_layers(kind)
    if layers:
        return _draw_layers(image, layers)

    renderer = orbital_renderer()
    spec = renderer.spec_for_kind(kind)
    painter = QPainter(image)
    anchor0, anchor1 = renderer.canonical_anchors(spec, QRectF(0.0, 0.0, float(CANVAS_SIZE), float(CANVAS_SIZE)))
    renderer.paint_glyph(painter, spec, anchor0, anchor1)
    painter.end()
    return image


def _render_current_crop(kind: str) -> QImage:
    image = QImage(CANVAS_SIZE, CANVAS_SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(BG)
    renderer = orbital_renderer()
    spec = renderer.spec_for_kind(kind)
    painter = QPainter(image)
    anchor0, anchor1 = renderer.canonical_anchors(spec, QRectF(0.0, 0.0, float(CANVAS_SIZE), float(CANVAS_SIZE)))
    renderer.paint_glyph(painter, spec, anchor0, anchor1)
    painter.end()
    return image


def _mask_from_image(image: QImage) -> set[tuple[int, int]]:
    mask: set[tuple[int, int]] = set()
    for y in range(image.height()):
        for x in range(image.width()):
            color = image.pixelColor(x, y)
            if color.alpha() <= 8:
                continue
            if color.red() > 245 and color.green() > 245 and color.blue() > 245:
                continue
            mask.add((x, y))
    return mask


def _bbox(mask: set[tuple[int, int]]) -> tuple[int, int, int, int]:
    xs = [x for x, _ in mask]
    ys = [y for _, y in mask]
    return min(xs), min(ys), max(xs), max(ys)


def _centroid(mask: set[tuple[int, int]]) -> tuple[float, float]:
    if not mask:
        return 0.0, 0.0
    sx = sum(x for x, _ in mask)
    sy = sum(y for _, y in mask)
    count = float(len(mask))
    return sx / count, sy / count


def _translate_mask(mask: set[tuple[int, int]], dx: int, dy: int) -> set[tuple[int, int]]:
    translated: set[tuple[int, int]] = set()
    for x, y in mask:
        nx = x + dx
        ny = y + dy
        if 0 <= nx < CANVAS_SIZE and 0 <= ny < CANVAS_SIZE:
            translated.add((nx, ny))
    return translated


def _align_current_mask(reference_mask: set[tuple[int, int]], current_mask: set[tuple[int, int]]) -> tuple[set[tuple[int, int]], tuple[float, float]]:
    rcx, rcy = _centroid(reference_mask)
    ccx, ccy = _centroid(current_mask)
    dx = int(round(rcx - ccx))
    dy = int(round(rcy - ccy))
    return _translate_mask(current_mask, dx, dy), (rcx - ccx, rcy - ccy)


def _overlay(reference_mask: set[tuple[int, int]], current_mask: set[tuple[int, int]]) -> QImage:
    image = QImage(CANVAS_SIZE, CANVAS_SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(BG)
    for y in range(CANVAS_SIZE):
        for x in range(CANVAS_SIZE):
            ref = (x, y) in reference_mask
            cur = (x, y) in current_mask
            if ref and cur:
                image.setPixelColor(x, y, QColor("#F6C544"))
            elif ref:
                image.setPixelColor(x, y, QColor("#D92B2B"))
            elif cur:
                image.setPixelColor(x, y, QColor("#2DA44E"))
    return image


def _compose_sheet(kind_to_image: dict[str, QImage], path: Path) -> QImage:
    image = QImage(ORBITAL_PALETTE_MODEL.image_size, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor("#F2F2F2"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for row in range(ORBITAL_PALETTE_MODEL.rows):
        for col in range(ORBITAL_PALETTE_MODEL.columns):
            cell = QRectF(
                col * ORBITAL_PALETTE_MODEL.cell_width,
                row * ORBITAL_PALETTE_MODEL.cell_height,
                ORBITAL_PALETTE_MODEL.cell_width,
                ORBITAL_PALETTE_MODEL.cell_height,
            )
            painter.fillRect(cell, QColor("#FFFFFF"))
            painter.setPen(QPen(QColor("#D8D8D8"), 1.0))
            painter.drawRect(cell.adjusted(0.0, 0.0, -1.0, -1.0))
            kind = ORBITAL_PALETTE_MODEL.kind_at(row, col)
            if not kind:
                continue
            crop = kind_to_image[kind]
            inner = cell.adjusted(
                ORBITAL_PALETTE_MODEL.inner_padding,
                ORBITAL_PALETTE_MODEL.inner_padding,
                -ORBITAL_PALETTE_MODEL.inner_padding,
                -ORBITAL_PALETTE_MODEL.inner_padding,
            )
            painter.drawImage(inner, crop)
    painter.end()
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(str(path))
    return image


def _compose_normalized_sheet(kind_to_image: dict[str, QImage], path: Path) -> QImage:
    image = QImage(
        ORBITAL_PALETTE_MODEL.columns * CANVAS_SIZE,
        ORBITAL_PALETTE_MODEL.rows * CANVAS_SIZE,
        QImage.Format.Format_ARGB32_Premultiplied,
    )
    image.fill(QColor("#F2F2F2"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for row in range(ORBITAL_PALETTE_MODEL.rows):
        for col in range(ORBITAL_PALETTE_MODEL.columns):
            cell = QRectF(
                float(col * CANVAS_SIZE),
                float(row * CANVAS_SIZE),
                float(CANVAS_SIZE),
                float(CANVAS_SIZE),
            )
            painter.fillRect(cell, QColor("#FFFFFF"))
            painter.setPen(QPen(QColor("#D8D8D8"), 1.0))
            painter.drawRect(cell.adjusted(0.0, 0.0, -1.0, -1.0))
            kind = ORBITAL_PALETTE_MODEL.kind_at(row, col)
            if not kind:
                continue
            painter.drawImage(cell.adjusted(1.0, 1.0, -1.0, -1.0), kind_to_image[kind])
    painter.end()
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(str(path))
    return image


def _cell_rect(kind: str, cell_size: int) -> QRectF:
    for row in range(ORBITAL_PALETTE_MODEL.rows):
        for col in range(ORBITAL_PALETTE_MODEL.columns):
            if ORBITAL_PALETTE_MODEL.kind_at(row, col) == kind:
                return QRectF(
                    float(col * cell_size),
                    float(row * cell_size),
                    float(cell_size),
                    float(cell_size),
                )
    raise KeyError(kind)


def _report_line(
    kind: str,
    iou: float,
    bbox_ratio_w: float,
    bbox_ratio_h: float,
    aspect_ratio_delta: float,
    centroid_px: float,
    extra_pct: float,
    missing_pct: float,
) -> str:
    return (
        f"{kind:<24} "
        f"IoU={iou:0.4f} "
        f"dW={bbox_ratio_w:0.3f} "
        f"dH={bbox_ratio_h:0.3f} "
        f"dAR={aspect_ratio_delta:0.3f} "
        f"dC={centroid_px:0.2f}px "
        f"extra={extra_pct:0.2f}% "
        f"missing={missing_pct:0.2f}%"
    )


def main() -> int:
    app = QApplication.instance() or QApplication([])
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_CROPS_DIR.mkdir(parents=True, exist_ok=True)
    CURRENT_CROPS_DIR.mkdir(parents=True, exist_ok=True)
    OVERLAY_DIR.mkdir(parents=True, exist_ok=True)

    reference_crops: dict[str, QImage] = {}
    current_crops: dict[str, QImage] = {}
    overlay_crops: dict[str, QImage] = {}
    report_lines: list[str] = []
    ious: list[float] = []

    for kind in ORBITAL_PALETTE_MODEL.entries:
        current = _render_current_crop(kind)
        current.save(str(CURRENT_CROPS_DIR / f"{kind}.png"))
        current_crops[kind] = current

    reference_source_crops = {kind: _render_reference_crop(kind) for kind in ORBITAL_PALETTE_MODEL.entries}
    for kind, reference in reference_source_crops.items():
        reference.save(str(REFERENCE_CROPS_DIR / f"{kind}.png"))
        reference_crops[kind] = QImage(str(REFERENCE_CROPS_DIR / f"{kind}.png")).convertToFormat(
            QImage.Format.Format_ARGB32
        )

    _compose_normalized_sheet(reference_crops, REFERENCE_SHEET)
    _compose_sheet(current_crops, GRID_PREVIEW_PATH)

    for kind in ORBITAL_PALETTE_MODEL.entries:
        reference_mask = _mask_from_image(reference_crops[kind])
        current_mask = _mask_from_image(current_crops[kind])
        aligned_current_mask, centroid_delta = _align_current_mask(reference_mask, current_mask)

        intersection = len(reference_mask & aligned_current_mask)
        union = len(reference_mask | aligned_current_mask) or 1
        iou = intersection / union
        ious.append(iou)

        ref_bbox = _bbox(reference_mask)
        cur_bbox = _bbox(aligned_current_mask)
        ref_w = max(1, ref_bbox[2] - ref_bbox[0] + 1)
        ref_h = max(1, ref_bbox[3] - ref_bbox[1] + 1)
        cur_w = max(1, cur_bbox[2] - cur_bbox[0] + 1)
        cur_h = max(1, cur_bbox[3] - cur_bbox[1] + 1)

        bbox_ratio_w = abs(cur_w - ref_w) / float(ref_w)
        bbox_ratio_h = abs(cur_h - ref_h) / float(ref_h)
        ref_aspect = ref_w / float(ref_h)
        cur_aspect = cur_w / float(cur_h)
        aspect_ratio_delta = abs(cur_aspect - ref_aspect) / float(ref_aspect or 1.0)
        centroid_px = math.hypot(*centroid_delta)
        extra_pct = (len(aligned_current_mask - reference_mask) / float(len(reference_mask) or 1)) * 100.0
        missing_pct = (len(reference_mask - aligned_current_mask) / float(len(reference_mask) or 1)) * 100.0

        overlay = _overlay(reference_mask, aligned_current_mask)
        overlay.save(str(OVERLAY_DIR / f"{kind}.png"))
        overlay_crops[kind] = overlay
        report_lines.append(
            _report_line(
                kind,
                iou,
                bbox_ratio_w,
                bbox_ratio_h,
                aspect_ratio_delta,
                centroid_px,
                extra_pct,
                missing_pct,
            )
        )

    summary = f"average IoU={sum(ious) / float(len(ious) or 1):0.4f}"
    _compose_normalized_sheet(overlay_crops, OVERLAY_SHEET_PATH)
    REPORT_PATH.write_text(
        "\n".join(
            [
                "reference=tests/data/orbitals/reference_sheet.png (redibujo vectorial fijo basado en la referencia del usuario)",
                summary,
                "",
            ]
            + report_lines
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary)
    print(f"report: {REPORT_PATH}")
    print(f"reference sheet: {REFERENCE_SHEET}")
    print(f"grid preview: {GRID_PREVIEW_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
