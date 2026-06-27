"""Verificación visual de la paleta orbital contra un baseline PNG."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtGui import QImage
from PyQt6.QtWidgets import QApplication

from chemuson.gui.orbitals import ORBITAL_PALETTE_MODEL, render_orbital_palette_image


_REPO_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_BASELINE = _REPO_ROOT / "assets" / "baseline" / "orbitals_palette.png"
_FOCAL_KINDS = {
    "sp3_shaded",
    "dz2_shaded",
    "sigma_bonding_shaded",
    "pi_bonding_shaded",
}


@dataclass(frozen=True)
class CellDiff:
    row: int
    column: int
    kind: str
    bad_pixels: int
    worst: int
    total_pixels: int

    @property
    def ratio(self) -> float:
        if self.total_pixels <= 0:
            return 0.0
        return self.bad_pixels / self.total_pixels


def _ensure_qapp() -> QApplication:
    app = QApplication.instance()
    return app if app is not None else QApplication(sys.argv)


def _load_image(path: Path) -> QImage:
    image = QImage(str(path)).convertToFormat(QImage.Format.Format_ARGB32)
    if image.isNull():
        raise FileNotFoundError(f"No se pudo cargar la imagen: {path}")
    return image


def _save_temp_render(image: QImage) -> Path:
    tmp = tempfile.NamedTemporaryFile(prefix="orbitals_palette_", suffix=".png", delete=False)
    tmp_path = Path(tmp.name)
    tmp.close()
    if not image.save(str(tmp_path)):
        raise RuntimeError(f"No se pudo guardar el render temporal: {tmp_path}")
    return tmp_path


def _pixel_delta(actual: QImage, expected: QImage, *, pixel_tolerance: int) -> tuple[int, int]:
    bad_pixels = 0
    worst = 0
    for y in range(actual.height()):
        for x in range(actual.width()):
            a = actual.pixelColor(x, y)
            b = expected.pixelColor(x, y)
            diff = max(
                abs(a.red() - b.red()),
                abs(a.green() - b.green()),
                abs(a.blue() - b.blue()),
                abs(a.alpha() - b.alpha()),
            )
            worst = max(worst, diff)
            if diff > pixel_tolerance:
                bad_pixels += 1
    return bad_pixels, worst


def _cell_delta(
    actual: QImage,
    expected: QImage,
    *,
    row: int,
    column: int,
    pixel_tolerance: int,
) -> tuple[int, int]:
    cell = ORBITAL_PALETTE_MODEL
    x0 = column * cell.cell_width
    y0 = row * cell.cell_height
    x1 = x0 + cell.cell_width
    y1 = y0 + cell.cell_height
    bad_pixels = 0
    worst = 0
    for y in range(y0, y1):
        for x in range(x0, x1):
            a = actual.pixelColor(x, y)
            b = expected.pixelColor(x, y)
            diff = max(
                abs(a.red() - b.red()),
                abs(a.green() - b.green()),
                abs(a.blue() - b.blue()),
                abs(a.alpha() - b.alpha()),
            )
            worst = max(worst, diff)
            if diff > pixel_tolerance:
                bad_pixels += 1
    return bad_pixels, worst


def _collect_cell_diffs(actual: QImage, expected: QImage, *, pixel_tolerance: int) -> list[CellDiff]:
    diffs: list[CellDiff] = []
    for row in range(ORBITAL_PALETTE_MODEL.rows):
        for column in range(ORBITAL_PALETTE_MODEL.columns):
            kind = ORBITAL_PALETTE_MODEL.kind_at(row, column)
            if not kind:
                continue
            bad_pixels, worst = _cell_delta(
                actual,
                expected,
                row=row,
                column=column,
                pixel_tolerance=pixel_tolerance,
            )
            diffs.append(
                CellDiff(
                    row=row + 1,
                    column=column + 1,
                    kind=kind,
                    bad_pixels=bad_pixels,
                    worst=worst,
                    total_pixels=ORBITAL_PALETTE_MODEL.cell_width * ORBITAL_PALETTE_MODEL.cell_height,
                )
            )
    return diffs


def _print_report(
    *,
    baseline_path: Path,
    temp_render_path: Path,
    overall_bad_pixels: int,
    overall_worst: int,
    overall_ratio: float,
    diffs: list[CellDiff],
) -> None:
    print(f"Baseline: {baseline_path}")
    print(f"Render temporal: {temp_render_path}")
    print(
        "Paleta completa: "
        f"bad_pixels={overall_bad_pixels} worst={overall_worst} ratio={overall_ratio:.4%}"
    )

    print("\nCasos focales:")
    for diff in diffs:
        if diff.kind not in _FOCAL_KINDS:
            continue
        print(
            f"  fila={diff.row} col={diff.column} kind={diff.kind:<22} "
            f"bad={diff.bad_pixels:4d} worst={diff.worst:3d} ratio={diff.ratio:.4%}"
        )

    changed = [diff for diff in diffs if diff.bad_pixels or diff.worst]
    if not changed:
        print("\nSin diferencias por celda.")
        return

    print("\nDiferencias por celda:")
    for diff in sorted(changed, key=lambda item: (-item.bad_pixels, -item.worst, item.row, item.column)):
        marker = "FOCAL" if diff.kind in _FOCAL_KINDS else "     "
        print(
            f"  {marker} fila={diff.row} col={diff.column} kind={diff.kind:<22} "
            f"bad={diff.bad_pixels:4d} worst={diff.worst:3d} ratio={diff.ratio:.4%}"
        )


def verify_orbitals_rendering(
    *,
    baseline_path: Path,
    pixel_tolerance: int = 6,
    max_worst: int = 24,
    max_bad_ratio: float = 0.004,
) -> int:
    _ensure_qapp()

    if not baseline_path.exists():
        print(f"Baseline faltante: {baseline_path}", file=sys.stderr)
        return 2

    actual = render_orbital_palette_image(ORBITAL_PALETTE_MODEL).convertToFormat(
        QImage.Format.Format_ARGB32
    )
    temp_render_path = _save_temp_render(actual)
    expected = _load_image(baseline_path)

    if actual.size() != expected.size():
        print(
            "Tamaño distinto entre render y baseline: "
            f"{actual.width()}x{actual.height()} vs {expected.width()}x{expected.height()}",
            file=sys.stderr,
        )
        print(f"Render temporal: {temp_render_path}", file=sys.stderr)
        return 2

    overall_bad_pixels, overall_worst = _pixel_delta(
        actual,
        expected,
        pixel_tolerance=pixel_tolerance,
    )
    total_pixels = actual.width() * actual.height()
    overall_ratio = overall_bad_pixels / total_pixels if total_pixels else 0.0
    diffs = _collect_cell_diffs(actual, expected, pixel_tolerance=pixel_tolerance)

    _print_report(
        baseline_path=baseline_path,
        temp_render_path=temp_render_path,
        overall_bad_pixels=overall_bad_pixels,
        overall_worst=overall_worst,
        overall_ratio=overall_ratio,
        diffs=diffs,
    )

    ok = overall_worst <= max_worst and overall_ratio <= max_bad_ratio
    print("\nRESULTADO:", "OK" if ok else "FAIL")
    return 0 if ok else 1


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=_DEFAULT_BASELINE,
        help=f"PNG baseline a comparar (default: {_DEFAULT_BASELINE})",
    )
    parser.add_argument(
        "--pixel-tolerance",
        type=int,
        default=6,
        help="Diferencia máxima por canal antes de contar un pixel como distinto.",
    )
    parser.add_argument(
        "--max-worst",
        type=int,
        default=24,
        help="Peor diferencia por pixel permitida.",
    )
    parser.add_argument(
        "--max-bad-ratio",
        type=float,
        default=0.004,
        help="Ratio máximo de pixeles fuera de tolerancia permitido.",
    )
    return parser


def main() -> int:
    parser = _build_arg_parser()
    args = parser.parse_args()
    return verify_orbitals_rendering(
        baseline_path=args.baseline,
        pixel_tolerance=args.pixel_tolerance,
        max_worst=args.max_worst,
        max_bad_ratio=args.max_bad_ratio,
    )


if __name__ == "__main__":
    raise SystemExit(main())
