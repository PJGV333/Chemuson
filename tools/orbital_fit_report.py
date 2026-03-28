"""Preview y reporte geométrico opcional para orbitales parametrizados."""

from __future__ import annotations

import argparse
import os
import statistics
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from PyQt6.QtCore import QRect, QRectF, Qt
from PyQt6.QtGui import QColor, QImage, QPainter, QPen
from PyQt6.QtWidgets import QApplication

from chemuson.gui.orbitals import (
    ORBITAL_PALETTE_MODEL,
    ORBITAL_PRESET_CONFIG_PATH,
    ORBITAL_SPECS,
    Dz2Params,
    FOrbitalParams,
    PiBondingParams,
    HybridOrbitalParams,
    OrbitalStyle,
    TorusParams,
    build_orbital_renderer,
    load_orbital_presets_payload,
)


REFERENCE_SHEET = ROOT / "tests" / "data" / "orbitals" / "reference_sheet.png"
OUTPUT_DIR = ROOT / "tests" / "data" / "orbitals" / "fit_report"
REFERENCE_CROPS_DIR = OUTPUT_DIR / "reference_crops"
CURRENT_CROPS_DIR = OUTPUT_DIR / "current_crops"
OVERLAY_DIR = OUTPUT_DIR / "overlay"
REPORT_PATH = OUTPUT_DIR / "report.txt"
OVERLAY_SHEET_PATH = OUTPUT_DIR / "overlay_sheet.png"
PALETTE_PREVIEW_PATH = ROOT / "tests" / "data" / "orbitals" / "palette_preview.png"

CANVAS_SIZE = 96
BG = QColor("#FFFFFF")
_SYMMETRY_AXES = {
    "s": ("vertical", "horizontal"),
    "p": ("vertical", "horizontal"),
    "d": ("vertical", "horizontal"),
    "dz2": ("vertical", "horizontal"),
    "torus": ("vertical", "horizontal"),
    "pi_bonding": ("vertical", "horizontal"),
    "fz3": ("vertical",),
}


def _ensure_qapp() -> QApplication:
    app = QApplication.instance()
    return app if app is not None else QApplication([])


def _load_image(path: Path) -> QImage:
    image = QImage(str(path)).convertToFormat(QImage.Format.Format_ARGB32)
    if image.isNull():
        raise FileNotFoundError(f"No se pudo cargar la imagen: {path}")
    return image


def _reference_cell_size(sheet: QImage) -> tuple[int, int]:
    return sheet.width() // ORBITAL_PALETTE_MODEL.columns, sheet.height() // ORBITAL_PALETTE_MODEL.rows


def _cell_position(kind: str) -> tuple[int, int]:
    for row in range(ORBITAL_PALETTE_MODEL.rows):
        for column in range(ORBITAL_PALETTE_MODEL.columns):
            if ORBITAL_PALETTE_MODEL.kind_at(row, column) == kind:
                return row, column
    raise KeyError(kind)


def _crop_reference(kind: str, sheet: QImage) -> QImage:
    cell_w, cell_h = _reference_cell_size(sheet)
    row, column = _cell_position(kind)
    crop = sheet.copy(QRect(column * cell_w + 1, row * cell_h + 1, cell_w - 2, cell_h - 2))
    return crop.scaled(
        CANVAS_SIZE,
        CANVAS_SIZE,
        Qt.AspectRatioMode.IgnoreAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    ).convertToFormat(QImage.Format.Format_ARGB32)


def _render_current_crop(renderer, kind: str) -> QImage:
    image = QImage(CANVAS_SIZE, CANVAS_SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(BG)
    spec = renderer.spec_for_kind(kind)
    painter = QPainter(image)
    anchor0, anchor1 = renderer.canonical_anchors(spec, QRectF(0.0, 0.0, float(CANVAS_SIZE), float(CANVAS_SIZE)))
    renderer.paint_glyph(painter, spec, anchor0, anchor1)
    painter.end()
    return image


def _render_silhouette(renderer, kind: str, style: OrbitalStyle) -> QImage:
    image = QImage(CANVAS_SIZE, CANVAS_SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(BG)
    base = ORBITAL_SPECS[kind]
    spec = type(base)(
        kind=base.kind,
        orbital_type=base.orbital_type,
        glyph_id=base.glyph_id,
        style=style,
        label=base.label,
        canvas_extent_scale=base.canvas_extent_scale,
        stroke_shaded_lobes=base.stroke_shaded_lobes,
    )
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    anchor0, anchor1 = renderer.canonical_anchors(spec, QRectF(0.0, 0.0, float(CANVAS_SIZE), float(CANVAS_SIZE)))
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor("#111111"))
    for layer in renderer.transformed_layers(spec, anchor0, anchor1):
        painter.drawPath(layer.path)
    painter.end()
    return image


def _mask_from_image(image: QImage) -> set[tuple[int, int]]:
    mask: set[tuple[int, int]] = set()
    for y in range(image.height()):
        for x in range(image.width()):
            color = image.pixelColor(x, y)
            if color.alpha() <= 8:
                continue
            if color.red() > 246 and color.green() > 246 and color.blue() > 246:
                continue
            mask.add((x, y))
    return mask


def _bbox(mask: set[tuple[int, int]]) -> tuple[int, int, int, int]:
    if not mask:
        return 0, 0, 0, 0
    xs = [x for x, _ in mask]
    ys = [y for _, y in mask]
    return min(xs), min(ys), max(xs), max(ys)


def _bbox_dims(mask: set[tuple[int, int]]) -> tuple[int, int]:
    left, top, right, bottom = _bbox(mask)
    return max(0, right - left + 1), max(0, bottom - top + 1)


def _centroid(mask: set[tuple[int, int]]) -> tuple[float, float]:
    if not mask:
        return 0.0, 0.0
    sx = sum(x for x, _ in mask)
    sy = sum(y for _, y in mask)
    count = float(len(mask))
    return sx / count, sy / count


def _occupancy(mask: set[tuple[int, int]]) -> float:
    return len(mask) / float(CANVAS_SIZE * CANVAS_SIZE)


def _iou(mask_a: set[tuple[int, int]], mask_b: set[tuple[int, int]]) -> float:
    union = len(mask_a | mask_b)
    if not union:
        return 1.0
    return len(mask_a & mask_b) / float(union)


def _symmetry_error(mask: set[tuple[int, int]], axis: str) -> float:
    if not mask:
        return 0.0
    mismatches = 0
    total = 0
    for x, y in mask:
        if axis == "vertical":
            mx, my = CANVAS_SIZE - 1 - x, y
        elif axis == "horizontal":
            mx, my = x, CANVAS_SIZE - 1 - y
        else:
            raise ValueError(axis)
        total += 1
        if (mx, my) not in mask:
            mismatches += 1
    return mismatches / float(total or 1)


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


def _silhouette_consistency(renderer, family: str) -> tuple[float, float, float]:
    outline = _mask_from_image(_render_silhouette(renderer, f"{family}_outline", OrbitalStyle.OUTLINE))
    shaded = _mask_from_image(_render_silhouette(renderer, f"{family}_shaded", OrbitalStyle.SHADED))
    solid = _mask_from_image(_render_silhouette(renderer, f"{family}_solid", OrbitalStyle.SOLID))
    return (
        _iou(outline, shaded),
        _iou(outline, solid),
        _iou(shaded, solid),
    )


def _format_optional(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.4f}"


def _report_line(kind: str, metrics: dict[str, str]) -> str:
    parts = [f"{kind:<24}"]
    for key, value in metrics.items():
        parts.append(f"{key}={value}")
    return " ".join(parts)


def _kind_family(kind: str) -> str:
    return kind.rsplit("_", 1)[0]


def _family_metric_strings(renderer, family: str) -> dict[str, str]:
    preset = renderer.family_preset(family)
    params = preset.params
    metrics: dict[str, str] = {}

    if isinstance(params, TorusParams):
        metrics["torus_ratio"] = f"{params.torus_inner_width_ratio:.3f}x{params.torus_inner_height_ratio:.3f}"
    elif isinstance(params, Dz2Params):
        metrics["torus_ratio"] = f"{params.torus_inner_width_ratio:.3f}x{params.torus_inner_height_ratio:.3f}"
        metrics["ring_offset"] = f"{params.torus_offset_y:.4f}"
    elif isinstance(params, PiBondingParams) and params.ring is not None:
        metrics["torus_ratio"] = (
            f"{params.ring.torus_inner_width_ratio:.3f}x"
            f"{params.ring.torus_inner_height_ratio:.3f}"
        )
    elif isinstance(params, HybridOrbitalParams):
        metrics["major_ar"] = f"{params.major_lobe.width / max(params.major_lobe.height, 1.0):.3f}"
        metrics["minor_ar"] = f"{params.minor_lobe.width / max(params.minor_lobe.height, 1.0):.3f}"
    elif isinstance(params, FOrbitalParams):
        metrics["lobes"] = str(len(params.lobes))
        if params.torus is not None:
            metrics["torus_ratio"] = (
                f"{params.torus.torus_inner_width_ratio:.3f}x"
                f"{params.torus.torus_inner_height_ratio:.3f}"
            )
    return metrics


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=ORBITAL_PRESET_CONFIG_PATH,
        help=f"Archivo de presets a inspeccionar (default: {ORBITAL_PRESET_CONFIG_PATH})",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=REFERENCE_SHEET if REFERENCE_SHEET.exists() else None,
        help="Hoja raster opcional para overlay. Si se omite, el reporte no fuerza comparacion externa.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    _ensure_qapp()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_CROPS_DIR.mkdir(parents=True, exist_ok=True)
    CURRENT_CROPS_DIR.mkdir(parents=True, exist_ok=True)
    OVERLAY_DIR.mkdir(parents=True, exist_ok=True)

    payload = load_orbital_presets_payload(args.config, include_docs=True)
    renderer = build_orbital_renderer(payload)

    reference_sheet = None
    if args.reference is not None:
        if not args.reference.exists():
            raise FileNotFoundError(f"Referencia no encontrada: {args.reference}")
        reference_sheet = _load_image(args.reference)

    reference_crops: dict[str, QImage] = {}
    current_crops: dict[str, QImage] = {}
    overlay_crops: dict[str, QImage] = {}
    report_lines: list[str] = []
    reference_ious: list[float] = []

    for kind in ORBITAL_PALETTE_MODEL.entries:
        current = _render_current_crop(renderer, kind)
        current.save(str(CURRENT_CROPS_DIR / f"{kind}.png"))
        current_crops[kind] = current

        if reference_sheet is not None:
            reference = _crop_reference(kind, reference_sheet)
            reference.save(str(REFERENCE_CROPS_DIR / f"{kind}.png"))
            reference_crops[kind] = reference

    palette = renderer.render_palette_image(ORBITAL_PALETTE_MODEL)
    palette.save(str(PALETTE_PREVIEW_PATH))

    for kind in ORBITAL_PALETTE_MODEL.entries:
        family = _kind_family(kind)
        current_mask = _mask_from_image(current_crops[kind])
        if family not in overlay_crops:
            pass

        ref_iou: float | None = None
        if reference_sheet is not None:
            reference_mask = _mask_from_image(reference_crops[kind])
            overlay = _overlay(reference_mask, current_mask)
            overlay.save(str(OVERLAY_DIR / f"{kind}.png"))
            overlay_crops[kind] = overlay
            ref_iou = _iou(reference_mask, current_mask)
            reference_ious.append(ref_iou)
        else:
            current_crops[kind].save(str(OVERLAY_DIR / f"{kind}.png"))
            overlay_crops[kind] = current_crops[kind]

        bbox = _bbox(current_mask)
        bbox_w, bbox_h = _bbox_dims(current_mask)
        centroid_x, centroid_y = _centroid(current_mask)
        centroid_dx = centroid_x - ((CANVAS_SIZE - 1) * 0.5)
        centroid_dy = centroid_y - ((CANVAS_SIZE - 1) * 0.5)
        occupancy = _occupancy(current_mask) * 100.0
        aspect_ratio = bbox_w / float(bbox_h or 1)

        sym_v: float | None = None
        sym_h: float | None = None
        for axis in _SYMMETRY_AXES.get(family, ()):
            error = _symmetry_error(
                _mask_from_image(_render_silhouette(renderer, f"{family}_outline", OrbitalStyle.OUTLINE)),
                axis,
            )
            if axis == "vertical":
                sym_v = error
            elif axis == "horizontal":
                sym_h = error

        sil_outline_shaded, sil_outline_solid, sil_shaded_solid = _silhouette_consistency(renderer, family)
        metrics: dict[str, str] = {
            "bbox": f"{bbox[0]},{bbox[1]},{bbox_w},{bbox_h}",
            "centroid": f"{centroid_dx:+.2f},{centroid_dy:+.2f}",
            "occupancy": f"{occupancy:.2f}%",
            "aspect": f"{aspect_ratio:.3f}",
            "ref_iou": _format_optional(ref_iou),
            "sym_v": _format_optional(sym_v),
            "sym_h": _format_optional(sym_h),
            "sil_o_s": f"{sil_outline_shaded:.4f}",
            "sil_o_f": f"{sil_outline_solid:.4f}",
            "sil_s_f": f"{sil_shaded_solid:.4f}",
        }
        metrics.update(_family_metric_strings(renderer, family))
        report_lines.append(_report_line(kind, metrics))

    _compose_normalized_sheet(overlay_crops, OVERLAY_SHEET_PATH)

    family_lines = []
    seen_families: list[str] = []
    for kind in ORBITAL_PALETTE_MODEL.entries:
        family = _kind_family(kind)
        if family in seen_families:
            continue
        seen_families.append(family)
        o_s, o_f, s_f = _silhouette_consistency(renderer, family)
        family_lines.append(
            f"{family:<18} outline/shaded={o_s:.4f} outline/solid={o_f:.4f} shaded/solid={s_f:.4f}"
        )

    average_ref_iou = statistics.fmean(reference_ious) if reference_ious else 0.0
    REPORT_PATH.write_text(
        "\n".join(
            [
                f"config={args.config}",
                f"reference_sheet={args.reference if reference_sheet is not None else 'none'}",
                "mode=preview_overlay_only",
                f"palette_preview={PALETTE_PREVIEW_PATH}",
                f"overlay_sheet={OVERLAY_SHEET_PATH}",
                f"average_ref_iou={average_ref_iou:.4f}",
                "",
                "Silhouette Consistency",
                *family_lines,
                "",
                "Per Glyph Metrics",
                *report_lines,
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"palette_preview={PALETTE_PREVIEW_PATH}")
    print(f"overlay_sheet={OVERLAY_SHEET_PATH}")
    print(f"report={REPORT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
