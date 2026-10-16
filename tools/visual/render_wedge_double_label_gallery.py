"""Render focused wedge/double-bond/label junctions through production Qt items.

Exploratory output defaults to /tmp. The gallery is deterministic, does not
move atoms after construction, and never changes chemical identity for testing.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QImage, QPainter, QPen
from PyQt6.QtWidgets import QApplication, QGraphicsRectItem, QGraphicsSimpleTextItem

from chemuson.core.model import BondStereo, BondStyle, MolGraph
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.templates import build_haworth_template

WIDTH = 1800
HEIGHT = 1760
COLUMNS = 3
CELL_WIDTH = 580
CELL_HEIGHT = 330
PAGE = QRectF(30.0, 34.0, WIDTH - 60.0, HEIGHT - 68.0)


CASES = (
    "hashed wedge + double, acute pi",
    "hashed wedge + double, near-parallel pi",
    "carbonyl O label · horizontal",
    "carbonyl O label · vertical",
    "carbonyl O label · diagonal",
    "solid wedge + F label",
    "hashed wedge + F label",
    "solid wedge + CH2OH label",
    "hashed wedge + CH2OH label",
    "Haworth hexagon · bold + solid wedge",
    "Haworth pentagon · bold + solid wedge",
    "control · ordinary double + O label",
    "Cl label · alkene horizontal",
    "N label · alkene vertical",
    "O label · alkene diagonal",
)


def _add_atom(
    canvas: ChemusonCanvas,
    element: str,
    x: float,
    y: float,
    *,
    explicit: bool = False,
    label_scale: float | None = None,
) -> int:
    return canvas.model.add_atom(
        element,
        x,
        y,
        is_explicit=explicit,
        label_scale=label_scale,
    ).id


def _add_bond(
    canvas: ChemusonCanvas,
    a1: int,
    a2: int,
    *,
    order: int = 1,
    style: BondStyle = BondStyle.PLAIN,
    stereo: BondStereo = BondStereo.NONE,
    pi_offset_sign: int | None = None,
) -> None:
    canvas.model.add_bond(
        a1,
        a2,
        order=order,
        style=style,
        stereo=stereo,
        pi_offset_sign=pi_offset_sign,
    )


def _add_hash_double(canvas: ChemusonCanvas, cx: float, cy: float, *, angle_deg: float, pi_sign: int) -> None:
    # Deliberately synthetic local junction: the renderer must keep sigma/pi
    # and the hashed mark profile visually distinct even in a crowded vertex.
    center = _add_atom(canvas, "C", cx, cy)
    angle = math.radians(180.0)
    tip = _add_atom(canvas, "C", cx + 96.0 * math.cos(angle), cy + 96.0 * math.sin(angle))
    _add_bond(canvas, tip, center, style=BondStyle.HASHED, stereo=BondStereo.DOWN)

    radians = math.radians(angle_deg)
    double_end = _add_atom(
        canvas,
        "C",
        cx + 96.0 * math.cos(radians),
        cy + 96.0 * math.sin(radians),
    )
    _add_bond(canvas, center, double_end, order=2, pi_offset_sign=pi_sign)
    continuation_angle = radians + math.radians(150.0)
    terminal = _add_atom(
        canvas,
        "C",
        canvas.model.get_atom(double_end).x + 96.0 * math.cos(continuation_angle),
        canvas.model.get_atom(double_end).y + 96.0 * math.sin(continuation_angle),
    )
    _add_bond(canvas, double_end, terminal)
    branch_angle = math.radians(80.0)
    branch = _add_atom(canvas, "C", cx + 96.0 * math.cos(branch_angle), cy + 96.0 * math.sin(branch_angle))
    _add_bond(canvas, center, branch)


def _add_carbonyl(canvas: ChemusonCanvas, cx: float, cy: float, angle_deg: float) -> None:
    oxygen = _add_atom(canvas, "O", cx, cy, explicit=True)
    radians = math.radians(angle_deg)
    carbon = _add_atom(canvas, "C", cx + 96.0 * math.cos(radians), cy + 96.0 * math.sin(radians))
    _add_bond(canvas, oxygen, carbon, order=2)
    # A second substituent gives the carbonyl end a directional neighbor while
    # leaving the O label as the measured endpoint.
    branch_angle = radians + math.radians(125.0)
    branch = _add_atom(
        canvas,
        "C",
        canvas.model.get_atom(carbon).x + 96.0 * math.cos(branch_angle),
        canvas.model.get_atom(carbon).y + 96.0 * math.sin(branch_angle),
    )
    _add_bond(canvas, carbon, branch)


def _add_alkene_label(canvas: ChemusonCanvas, cx: float, cy: float, label: str, angle_deg: float) -> None:
    alkene = _add_atom(canvas, "C", cx, cy)
    angle = math.radians(angle_deg)
    alkene_end = _add_atom(canvas, "C", cx + 96.0 * math.cos(angle), cy + 96.0 * math.sin(angle))
    _add_bond(canvas, alkene, alkene_end, order=2)
    continuation_angle = angle + math.radians(145.0)
    continuation = _add_atom(
        canvas,
        "C",
        canvas.model.get_atom(alkene_end).x + 96.0 * math.cos(continuation_angle),
        canvas.model.get_atom(alkene_end).y + 96.0 * math.sin(continuation_angle),
    )
    _add_bond(canvas, alkene_end, continuation)
    branch_angle = angle + math.radians(135.0)
    substituent = _add_atom(
        canvas,
        label,
        cx + 96.0 * math.cos(branch_angle),
        cy + 96.0 * math.sin(branch_angle),
        explicit=True,
    )
    _add_bond(canvas, alkene, substituent)


def _add_labeled_wedge(
    canvas: ChemusonCanvas,
    cx: float,
    cy: float,
    label: str,
    style: BondStyle,
    *,
    label_scale: float = 1.0,
) -> None:
    label_atom = _add_atom(canvas, label, cx, cy, explicit=True, label_scale=label_scale)
    other = _add_atom(canvas, "C", cx + 96.0, cy)
    stereo = BondStereo.UP if style == BondStyle.WEDGE else BondStereo.DOWN
    _add_bond(canvas, label_atom, other, style=style, stereo=stereo)


def _add_haworth_graph(canvas: ChemusonCanvas, graph: MolGraph, cx: float, cy: float) -> None:
    atom_map: dict[int, int] = {}
    xs = [atom.x for atom in graph.atoms.values()]
    ys = [atom.y for atom in graph.atoms.values()]
    graph_cx = (min(xs) + max(xs)) * 0.5
    graph_cy = (min(ys) + max(ys)) * 0.5
    for atom in graph.atoms.values():
        atom_map[atom.id] = canvas.model.add_atom(
            atom.element,
            cx + atom.x - graph_cx,
            cy + atom.y - graph_cy,
            charge=atom.charge,
            explicit_h=atom.explicit_h,
            is_explicit=atom.is_explicit,
            label_scale=atom.label_scale,
        ).id
    for bond in graph.bonds.values():
        _add_bond(
            canvas,
            atom_map[bond.a1_id],
            atom_map[bond.a2_id],
            order=bond.order,
            style=bond.style,
            stereo=bond.stereo,
        )


def _build_haworth_pentagon(bond_length: float = 78.0) -> MolGraph:
    graph = MolGraph()
    coords = (
        (0.0, -0.72 * bond_length),
        (0.72 * bond_length, -0.20 * bond_length),
        (0.45 * bond_length, 0.62 * bond_length),
        (-0.45 * bond_length, 0.62 * bond_length),
        (-0.72 * bond_length, -0.20 * bond_length),
    )
    elements = ("O", "C", "C", "C", "C")
    ring = [
        graph.add_atom(element, x, y, is_explicit=element != "C")
        for element, (x, y) in zip(elements, coords)
    ]
    for index in range(len(ring)):
        style = BondStyle.BOLD if index in {2, 3} else BondStyle.PLAIN
        graph.add_bond(ring[index].id, ring[(index + 1) % len(ring)].id, style=style)
    # Anomeric OH with a solid wedge, adjacent to the bold front edge.
    oh = graph.add_atom("OH", ring[3].x - 0.72 * bond_length, ring[3].y - 0.62 * bond_length, is_explicit=True)
    # Orient the wedge from substituent tip toward the ring so the broad end
    # meets both bold front edges at the lower-left ring atom.
    graph.add_bond(oh.id, ring[3].id, style=BondStyle.WEDGE, stereo=BondStereo.UP)
    return graph


def _add_case(canvas: ChemusonCanvas, title: str, cx: float, cy: float) -> None:
    if title == CASES[0]:
        _add_hash_double(canvas, cx, cy, angle_deg=35.0, pi_sign=1)
    elif title == CASES[1]:
        _add_hash_double(canvas, cx, cy, angle_deg=5.0, pi_sign=-1)
    elif title == CASES[2]:
        _add_carbonyl(canvas, cx, cy, 0.0)
    elif title == CASES[3]:
        _add_carbonyl(canvas, cx, cy, 90.0)
    elif title == CASES[4]:
        _add_carbonyl(canvas, cx, cy, 45.0)
    elif title == CASES[5]:
        _add_labeled_wedge(canvas, cx, cy, "F", BondStyle.WEDGE)
    elif title == CASES[6]:
        _add_labeled_wedge(canvas, cx, cy, "F", BondStyle.HASHED)
    elif title == CASES[7]:
        _add_labeled_wedge(canvas, cx, cy, "CH2OH", BondStyle.WEDGE, label_scale=1.15)
    elif title == CASES[8]:
        _add_labeled_wedge(canvas, cx, cy, "CH2OH", BondStyle.HASHED, label_scale=1.15)
    elif title == CASES[9]:
        graph = build_haworth_template(70.0, anomeric_up=True, bold_front=True)
        for bond in graph.bonds.values():
            atom_a = graph.get_atom(bond.a1_id)
            atom_b = graph.get_atom(bond.a2_id)
            if bond.a1_id == 4 and atom_b.element == "O" and atom_b.x < atom_a.x:
                # C3-OH belongs to a bold-front vertex. Keep the graph edge,
                # but point the wedge tip at the substituent and its broad end
                # into the ring junction.
                bond.a1_id, bond.a2_id = atom_b.id, atom_a.id
                bond.style = BondStyle.WEDGE
                bond.stereo = BondStereo.UP
                break
        _add_haworth_graph(canvas, graph, cx, cy)
    elif title == CASES[10]:
        _add_haworth_graph(canvas, _build_haworth_pentagon(), cx, cy)
    elif title == CASES[11]:
        _add_carbonyl(canvas, cx, cy, 180.0)
    elif title == CASES[12]:
        _add_alkene_label(canvas, cx, cy, "Cl", 0.0)
    elif title == CASES[13]:
        _add_alkene_label(canvas, cx, cy, "N", 90.0)
    elif title == CASES[14]:
        _add_alkene_label(canvas, cx, cy, "O", 45.0)
    else:  # pragma: no cover - CASES is exhaustive
        raise ValueError(f"Unknown gallery case: {title}")


def render_gallery(theme: str, output: Path) -> QImage:
    if theme not in {"light", "dark"}:
        raise ValueError(f"Unsupported theme: {theme}")
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

    heading = QGraphicsSimpleTextItem(f"ChemUSON · wedge / double / labels · {theme.upper()}")
    heading.setBrush(QBrush(QColor("#0F172A" if theme == "light" else "#F8FAFC")))
    heading.setFont(QFont("DejaVu Sans", 17, QFont.Weight.DemiBold))
    heading.setPos(42.0, 5.0)
    heading.setZValue(100.0)
    canvas.scene.addItem(heading)

    for index, title in enumerate(CASES):
        row, column = divmod(index, COLUMNS)
        x0 = PAGE.x() + column * CELL_WIDTH
        y0 = PAGE.y() + 42.0 + row * CELL_HEIGHT
        caption = QGraphicsSimpleTextItem(f"{index + 1:02d} · {title}")
        caption.setBrush(QBrush(QColor("#334155")))
        caption.setFont(QFont("DejaVu Sans", 10, QFont.Weight.Medium))
        caption.setPos(x0 + 8.0, y0 + 4.0)
        caption.setZValue(100.0)
        canvas.scene.addItem(caption)
        _add_case(canvas, title, x0 + CELL_WIDTH * 0.53, y0 + CELL_HEIGHT * 0.61)

    canvas._rebuild_items_from_model()
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


def compose_comparison(before_dir: Path, after_dir: Path, output: Path) -> None:
    """Compose equal-scale light/dark before/after panels for owner review."""
    app = QApplication.instance() or QApplication(sys.argv)
    image_height = HEIGHT + 44
    sheet = QImage(2 * WIDTH, 2 * image_height, QImage.Format.Format_ARGB32_Premultiplied)
    sheet.fill(QColor("#E2E8F0"))
    painter = QPainter(sheet)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.setFont(QFont("DejaVu Sans", 14, QFont.Weight.DemiBold))
    for row, theme in enumerate(("light", "dark")):
        for column, (phase, directory) in enumerate((("BEFORE", before_dir), ("AFTER", after_dir))):
            image_path = directory / f"wedge-double-labels-{theme}.png"
            image = QImage(str(image_path))
            if image.isNull() or image.size().width() != WIDTH or image.size().height() != HEIGHT:
                painter.end()
                raise ValueError(f"Missing or incorrectly sized gallery panel: {image_path}")
            x = column * WIDTH
            y = row * image_height
            header_color = QColor("#F1F5F9" if theme == "light" else "#1B263A")
            painter.fillRect(QRectF(x, y, WIDTH, 44), header_color)
            painter.setPen(QColor("#0F172A" if theme == "light" else "#F8FAFC"))
            painter.drawText(x + 18, y + 29, f"{theme.upper()} · {phase}")
            painter.drawImage(x, y + 44, image)
    painter.end()
    output.parent.mkdir(parents=True, exist_ok=True)
    if not sheet.save(str(output)):
        raise RuntimeError(f"Could not save comparison sheet: {output}")
    app.processEvents()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="both")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/chemuson-wedge-double-label-gallery"),
    )
    parser.add_argument("--compare-before", type=Path)
    parser.add_argument("--compare-after", type=Path)
    parser.add_argument("--comparison-output", type=Path)
    args = parser.parse_args(argv)
    themes = ("light", "dark") if args.theme == "both" else (args.theme,)
    for theme in themes:
        path = args.output_dir / f"wedge-double-labels-{theme}.png"
        render_gallery(theme, path)
        print(path)
    if any((args.compare_before, args.compare_after, args.comparison_output)):
        if not all((args.compare_before, args.compare_after, args.comparison_output)):
            parser.error("--compare-before, --compare-after, and --comparison-output must be set together")
        compose_comparison(args.compare_before, args.compare_after, args.comparison_output)
        print(args.comparison_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
