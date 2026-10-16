"""Render a deterministic Qt gallery of molecular bond junctions.

The gallery uses the production ChemusonCanvas/BondItem pipeline. It writes
exploratory PNGs to a caller-selected directory (default: /tmp), never changes
molecular coordinates after construction, and performs no chemical layout.
"""
from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QImage, QPainter, QPen
from PyQt6.QtWidgets import QApplication, QGraphicsRectItem, QGraphicsSimpleTextItem

from chemuson.core.model import BondStyle, BondStereo
from chemuson.gui.canvas import ChemusonCanvas


WIDTH = 1680
HEIGHT = 1160
COLUMNS = 4
ROWS = 5
CELL_WIDTH = 400
CELL_HEIGHT = 195
PAGE = QRectF(36.0, 36.0, WIDTH - 72.0, HEIGHT - 72.0)


@dataclass(frozen=True)
class Branch:
    angle_deg: float
    style: BondStyle = BondStyle.PLAIN
    order: int = 1
    stroke_px: float = 2.6
    element: str = "C"


@dataclass(frozen=True)
class GalleryCase:
    title: str
    branches: tuple[Branch, ...] = ()
    center_element: str = "C"
    charge: int = 0
    isotope: int | None = None
    stereo_style: BondStyle = BondStyle.WEDGE
    stereo_angle_deg: float = 180.0
    stereo_length: float = 50.0
    stereo_stroke: float = 2.6
    ring: str | None = None
    explicit_h: bool = False


def gallery_cases() -> tuple[GalleryCase, ...]:
    """Return fixed synthetic depictions; angles are measured at the junction."""
    return (
        GalleryCase("01 · simple + solid wedge · 30°", (Branch(30),)),
        GalleryCase(
            "02 · simple + hashed wedge · 60°",
            (Branch(60),),
            stereo_style=BondStyle.HASHED,
        ),
        GalleryCase("03 · sigma/pi double + solid wedge · 90°", (Branch(90, order=2),)),
        GalleryCase("04 · bold + solid wedge · 120°", (Branch(120, BondStyle.BOLD),)),
        GalleryCase(
            "05 · bold + hashed wedge · 150°",
            (Branch(150, BondStyle.BOLD),),
            stereo_style=BondStyle.HASHED,
        ),
        GalleryCase(
            "06 · sigma/pi double + hashed wedge · 30°",
            (Branch(30, order=2),),
            stereo_style=BondStyle.HASHED,
        ),
        GalleryCase("07 · solid wedge + hashed wedge", (Branch(0, BondStyle.HASHED),)),
        GalleryCase("08 · solid wedge + solid wedge", (Branch(0, BondStyle.WEDGE),)),
        GalleryCase(
            "09 · three bonds at a junction",
            (Branch(-35), Branch(35, BondStyle.BOLD)),
        ),
        GalleryCase(
            "10 · four bonds at a junction",
            (Branch(-65), Branch(0, BondStyle.BOLD), Branch(65)),
        ),
        GalleryCase("11 · near-collinear 178°", (Branch(178),)),
        GalleryCase(
            "12 · 30° fork · thin + thick strokes",
            (Branch(-15, stroke_px=1.3), Branch(15, stroke_px=4.0)),
            stereo_stroke=4.0,
        ),
        GalleryCase(
            "13 · charged N + explicit H",
            (Branch(35), Branch(130, BondStyle.BOLD)),
            center_element="N",
            charge=1,
            explicit_h=True,
        ),
        GalleryCase(
            "14 · O label at wedge junction",
            (Branch(60),),
            center_element="O",
        ),
        GalleryCase(
            "15 · S label + triple bond",
            (Branch(70, order=3),),
            center_element="S",
        ),
        GalleryCase(
            "16 · short scale · fine strokes",
            (Branch(120, stroke_px=1.3),),
            stereo_length=28.0,
            stereo_stroke=1.3,
        ),
        GalleryCase(
            "17 · long scale · bold strokes",
            (Branch(30, BondStyle.BOLD, stroke_px=4.0),),
            stereo_length=64.0,
            stereo_stroke=4.0,
        ),
        GalleryCase(
            "18 · charged N, label clearance",
            (Branch(150),),
            center_element="N",
            charge=-1,
            stereo_stroke=3.5,
        ),
        GalleryCase("19 · aromatic ring + solid wedge", ring="aromatic"),
        GalleryCase("20 · fused aromatic rings + wedge", ring="fused"),
    )


def _add_bond(canvas: ChemusonCanvas, a1: int, a2: int, branch: Branch) -> None:
    canvas.model.add_bond(
        a1,
        a2,
        order=branch.order,
        style=branch.style,
        stereo=(BondStereo.UP if branch.style == BondStyle.WEDGE else BondStereo.DOWN)
        if branch.style in {BondStyle.WEDGE, BondStyle.HASHED}
        else BondStereo.NONE,
        stroke_px=branch.stroke_px,
    )


def _add_aromatic_ring(canvas: ChemusonCanvas, cx: float, cy: float, radius: float, ring_id: int):
    coords = [
        (
            cx + radius * math.cos(math.radians(60.0 * i)),
            cy + radius * math.sin(math.radians(60.0 * i)),
        )
        for i in range(6)
    ]
    atoms = [canvas.model.add_atom("C", x, y) for x, y in coords]
    for index in range(6):
        a1 = atoms[index]
        a2 = atoms[(index + 1) % 6]
        order = 2 if index % 2 == 0 else 1
        canvas.model.add_bond(
            a1.id,
            a2.id,
            order=order,
            is_aromatic=True,
            display_order=order,
            ring_id=ring_id,
        )
    return atoms


def _add_case(canvas: ChemusonCanvas, case: GalleryCase, x0: float, y0: float, index: int) -> None:
    cx = x0 + CELL_WIDTH * 0.54
    cy = y0 + CELL_HEIGHT * 0.60
    if case.ring:
        radius = 35.0
        if case.ring == "aromatic":
            ring_atoms = _add_aromatic_ring(canvas, cx, cy, radius, index)
            root = ring_atoms[0]
            tip = canvas.model.add_atom("C", root.x - 39.0, root.y - 25.0)
            canvas.model.add_bond(
                tip.id,
                root.id,
                style=BondStyle.WEDGE,
                stereo=BondStereo.UP,
                stroke_px=case.stereo_stroke,
            )
            return

        # Two regular hexagons share one complete edge (a fused system).
        left = _add_aromatic_ring(canvas, cx - radius * 0.5, cy, radius, index)
        right = _add_aromatic_ring(canvas, cx + radius * 0.5, cy, radius, index + 100)
        # Merge the coincident shared-edge atoms into one molecular graph atom.
        # Vertices at +60/-60 degrees on the left ring coincide with
        # +120/+240 degrees on the right ring when the centers are one radius apart.
        shared_pairs = ((1, 2), (5, 4))
        for left_i, right_i in shared_pairs:
            keep, discard = left[left_i], right[right_i]
            for bond in list(canvas.model.bonds.values()):
                if bond.a1_id == discard.id:
                    bond.a1_id = keep.id
                if bond.a2_id == discard.id:
                    bond.a2_id = keep.id
            canvas.model.atoms.pop(discard.id, None)
        # The shared edge was inserted twice; retain one copy per endpoint pair.
        seen: set[frozenset[int]] = set()
        for bond_id, bond in list(canvas.model.bonds.items()):
            key = frozenset((bond.a1_id, bond.a2_id))
            if key in seen:
                canvas.model.bonds.pop(bond_id, None)
            else:
                seen.add(key)
        root = left[0]
        tip = canvas.model.add_atom("C", root.x - 40.0, root.y - 24.0)
        canvas.model.add_bond(tip.id, root.id, style=BondStyle.WEDGE, stereo=BondStereo.UP)
        return

    center = canvas.model.add_atom(
        case.center_element,
        cx,
        cy,
        formal_charge=case.charge,
        isotope=case.isotope,
        is_explicit=case.center_element != "C",
    )
    theta = math.radians(case.stereo_angle_deg)
    tip = canvas.model.add_atom(
        "C",
        cx - math.cos(theta) * case.stereo_length,
        cy - math.sin(theta) * case.stereo_length,
    )
    canvas.model.add_bond(
        tip.id,
        center.id,
        style=case.stereo_style,
        stereo=BondStereo.UP if case.stereo_style == BondStyle.WEDGE else BondStereo.DOWN,
        stroke_px=case.stereo_stroke,
    )

    for branch in case.branches:
        angle = math.radians(branch.angle_deg)
        end = canvas.model.add_atom(
            branch.element,
            cx + math.cos(angle) * case.stereo_length,
            cy + math.sin(angle) * case.stereo_length,
            is_explicit=branch.element != "C",
        )
        _add_bond(canvas, center.id, end.id, branch)
    if case.explicit_h:
        angle = math.radians(225.0)
        hydrogen = canvas.model.add_atom(
            "H",
            cx + math.cos(angle) * case.stereo_length,
            cy + math.sin(angle) * case.stereo_length,
            is_explicit=True,
        )
        canvas.model.add_bond(center.id, hydrogen.id)


def render_gallery(theme: str, output: Path, *, numbering: bool = False) -> QImage:
    """Render one light/dark context using the production Qt graphics items."""
    app = QApplication.instance() or QApplication(sys.argv)
    canvas = ChemusonCanvas()
    canvas.paper.setVisible(False)
    for grid_item in (canvas._grid_minor_item, canvas._grid_major_item):
        if grid_item is not None:
            grid_item.setVisible(False)

    background = QColor("#F1F5F9" if theme == "light" else "#1B263A")
    canvas.scene.setBackgroundBrush(QBrush(background))
    page = QGraphicsRectItem(PAGE)
    page.setBrush(QBrush(QColor("#FFFFFF")))
    page.setPen(QPen(QColor("#CBD5E1"), 1.0))
    page.setZValue(-20.0)
    canvas.scene.addItem(page)

    # Keep the heading on the theme surface; all molecular examples are on the
    # same white sheet in both themes, as in the application's page canvas.
    ink = QColor("#0F172A" if theme == "light" else "#F8FAFC")
    heading = QGraphicsSimpleTextItem(f"ChemUSON · bond junctions · {theme.upper()}")
    heading.setBrush(QBrush(ink))
    heading.setFont(QFont("DejaVu Sans", 15, QFont.Weight.DemiBold))
    heading.setPos(46.0, 7.0)
    heading.setZValue(100.0)
    canvas.scene.addItem(heading)

    for index, case in enumerate(gallery_cases()):
        row, column = divmod(index, COLUMNS)
        x0 = PAGE.x() + column * CELL_WIDTH
        y0 = PAGE.y() + 54.0 + row * CELL_HEIGHT
        caption = QGraphicsSimpleTextItem(case.title)
        caption.setBrush(QBrush(QColor("#334155")))
        caption.setFont(QFont("DejaVu Sans", 8))
        caption.setPos(x0 + 6.0, y0 + 2.0)
        caption.setZValue(100.0)
        canvas.scene.addItem(caption)
        _add_case(canvas, case, x0, y0, index + 1)

    canvas._rebuild_items_from_model()
    if numbering:
        canvas.set_numbering_enabled(True)
    canvas.scene.setSceneRect(0.0, 0.0, WIDTH, HEIGHT)
    image = QImage(WIDTH, HEIGHT, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(background)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    canvas.scene.render(
        painter,
        QRectF(0.0, 0.0, WIDTH, HEIGHT),
        QRectF(0.0, 0.0, WIDTH, HEIGHT),
        Qt.AspectRatioMode.IgnoreAspectRatio,
    )
    painter.end()
    output.parent.mkdir(parents=True, exist_ok=True)
    if not image.save(str(output)):
        raise RuntimeError(f"Could not save gallery image: {output}")
    canvas.close()
    app.processEvents()
    return image


def compose_comparison(
    before_dir: Path,
    after_dir: Path,
    output: Path,
    *,
    numbering: bool = False,
) -> None:
    """Compose one two-theme before/after reference sheet with Qt only."""
    app = QApplication.instance() or QApplication(sys.argv)
    panel_header = 44
    panel_height = HEIGHT + panel_header
    image = QImage(WIDTH * 2, panel_height * 2, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor("#E2E8F0"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.setFont(QFont("DejaVu Sans", 15, QFont.Weight.DemiBold))
    for row, theme in enumerate(("light", "dark")):
        for column, (stage, directory) in enumerate((("BEFORE", before_dir), ("AFTER", after_dir))):
            suffix = "-numbered" if numbering else ""
            path = directory / f"bond-junctions-{theme}{suffix}.png"
            panel = QImage(str(path))
            if panel.isNull() or panel.width() != WIDTH or panel.height() != HEIGHT:
                painter.end()
                raise ValueError(f"Missing or wrong-sized gallery panel: {path}")
            x = column * WIDTH
            y = row * panel_height
            header_color = QColor("#F8FAFC" if theme == "light" else "#1B263A")
            text_color = QColor("#0F172A" if theme == "light" else "#F8FAFC")
            painter.fillRect(x, y, WIDTH, panel_header, header_color)
            painter.setPen(text_color)
            painter.drawText(x + 18, y + 29, f"{stage} · {theme.upper()}")
            painter.drawImage(x, y + panel_header, panel)
    painter.end()
    output.parent.mkdir(parents=True, exist_ok=True)
    if not image.save(str(output)):
        raise RuntimeError(f"Could not save comparison sheet: {output}")
    app.processEvents()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="both")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/chemuson-bond-junction-gallery"),
        help="Output directory (defaults to /tmp; no files are written to the repository).",
    )
    parser.add_argument("--numbering", action="store_true", help="Include real atom-numbering overlays.")
    parser.add_argument("--before-dir", type=Path, help="Compose with matching before images from this directory.")
    parser.add_argument("--comparison-output", type=Path, help="Write one 2×2 before/after light/dark sheet.")
    args = parser.parse_args(argv)
    if bool(args.before_dir) != bool(args.comparison_output):
        parser.error("--before-dir and --comparison-output must be supplied together")
    themes = ("light", "dark") if args.theme == "both" else (args.theme,)
    suffix = "-numbered" if args.numbering else ""
    for theme in themes:
        path = args.output_dir / f"bond-junctions-{theme}{suffix}.png"
        render_gallery(theme, path, numbering=args.numbering)
        print(path)
    if args.before_dir and args.comparison_output:
        compose_comparison(
            args.before_dir,
            args.output_dir,
            args.comparison_output,
            numbering=args.numbering,
        )
        print(args.comparison_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
