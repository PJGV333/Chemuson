from __future__ import annotations

import hashlib
import math

import pytest
from PyQt6.QtGui import QColor, QImage
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import Atom, Bond, BondStyle, BondStereo
from chemuson.gui.bond_junctions import hashed_wedge_segments
from chemuson.gui.items import BondItem
from tools.visual.render_bond_junction_gallery import compose_comparison, render_gallery


@pytest.mark.parametrize("angle_deg", (30.0, 60.0, 90.0, 120.0, 150.0, 179.0))
@pytest.mark.parametrize("length,scale", ((24.0, 0.5), (50.0, 1.0), (72.0, 1.5)))
def test_hash_strokes_are_uniformly_spaced_and_grow_with_scale(angle_deg, length, scale):
    angle = math.radians(angle_deg)
    start = (17.0, -11.0)
    end = (start[0] + length * math.cos(angle), start[1] + length * math.sin(angle))
    count = 6
    min_width, max_width = 2.0 * scale, 10.0 * scale

    segments = hashed_wedge_segments(start, end, count, min_width, max_width)

    assert len(segments) == count
    ux, uy = math.cos(angle), math.sin(angle)
    centers = []
    widths = []
    for index, (p0, p1) in enumerate(segments):
        cx, cy = (p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5
        fraction = ((cx - start[0]) * ux + (cy - start[1]) * uy) / length
        centers.append(fraction)
        widths.append(math.hypot(p1[0] - p0[0], p1[1] - p0[1]))
        vx, vy = p1[0] - p0[0], p1[1] - p0[1]
        assert abs(vx * ux + vy * uy) < 1e-9
        assert fraction == pytest.approx((index + 0.5) / count, abs=1e-9)
        assert widths[-1] == pytest.approx(
            min_width + (max_width - min_width) * fraction,
            abs=1e-8,
        )

    expected_spacing = 1.0 / count
    assert [b - a for a, b in zip(centers, centers[1:])] == pytest.approx(
        [expected_spacing] * (count - 1), abs=1e-9
    )
    assert widths == sorted(widths)
    assert centers[0] == pytest.approx(0.5 / count)
    assert centers[-1] == pytest.approx(1.0 - 0.5 / count)


def test_hash_strokes_are_rotation_translation_and_repeat_deterministic():
    start = (-4.5, 17.25)
    end = (42.0, -5.75)
    expected = hashed_wedge_segments(start, end, 7, 1.2, 8.4)
    assert hashed_wedge_segments(start, end, 7, 1.2, 8.4) == expected

    angle = math.radians(73.0)
    c, s = math.cos(angle), math.sin(angle)
    rotated = tuple((x * c - y * s + 100.0, x * s + y * c - 80.0) for x, y in (start, end))
    actual_rotated = hashed_wedge_segments(rotated[0], rotated[1], 7, 1.2, 8.4)
    expected_rotated = tuple(
        (
            (p0[0] * c - p0[1] * s + 100.0, p0[0] * s + p0[1] * c - 80.0),
            (p1[0] * c - p1[1] * s + 100.0, p1[0] * s + p1[1] * c - 80.0),
        )
        for p0, p1 in expected
    )
    for actual_pair, expected_pair in zip(actual_rotated, expected_rotated):
        for actual, target in zip(actual_pair, expected_pair):
            assert actual == pytest.approx(target, abs=1e-9)


@pytest.mark.parametrize(
    "args",
    (
        ((0.0, 0.0), (0.0, 0.0), 6, 2.0, 10.0),
        ((0.0, 0.0), (10.0, 0.0), 0, 2.0, 10.0),
        ((0.0, 0.0), (10.0, 0.0), 6, 10.0, 2.0),
        ((0.0, 0.0), (10.0, 0.0), 6, 2.0, math.inf),
    ),
)
def test_hash_stroke_geometry_fails_closed_for_degenerate_inputs(args):
    assert hashed_wedge_segments(*args) == ()


def test_bond_item_uses_centered_hash_profile_without_mutating_graph():
    app = QApplication.instance() or QApplication([])
    atom1 = Atom(id=1, element="C", x=5.0, y=9.0)
    atom2 = Atom(id=2, element="C", x=65.0, y=9.0)
    bond = Bond(
        id=1,
        a1_id=atom1.id,
        a2_id=atom2.id,
        style=BondStyle.HASHED,
        stereo=BondStereo.DOWN,
        stroke_px=2.6,
    )
    snapshot = (
        (atom1.id, atom1.element, atom1.x, atom1.y),
        (atom2.id, atom2.element, atom2.x, atom2.y),
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo, bond.stroke_px),
    )

    item = BondItem(bond, atom1, atom2)

    assert item.path().elementCount() == 12
    centers = []
    for index in range(6):
        start = item.path().elementAt(index * 2)
        end = item.path().elementAt(index * 2 + 1)
        centers.append((start.x + end.x) * 0.5)
    assert centers == pytest.approx([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    assert snapshot == (
        (atom1.id, atom1.element, atom1.x, atom1.y),
        (atom2.id, atom2.element, atom2.x, atom2.y),
        (bond.id, bond.a1_id, bond.a2_id, bond.order, bond.style, bond.stereo, bond.stroke_px),
    )
    item.update_positions(atom1, atom2)
    assert item.path().elementCount() == 12
    app.processEvents()


def _orientation(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _proper_segments_intersect(a, b, c, d):
    epsilon = 1e-8
    return (
        _orientation(a, b, c) * _orientation(a, b, d) < -epsilon
        and _orientation(c, d, a) * _orientation(c, d, b) < -epsilon
    )


def test_solid_wedge_join_is_finite_and_non_self_intersecting_across_angles():
    app = QApplication.instance() or QApplication([])
    for angle_deg in (30.0, 60.0, 90.0, 120.0, 150.0, 179.0):
        for length, stroke in ((28.0, 1.3), (50.0, 2.6), (72.0, 4.0)):
            atom1 = Atom(id=1, element="C", x=-length, y=0.0)
            atom2 = Atom(id=2, element="C", x=0.0, y=0.0)
            bond = Bond(
                id=1,
                a1_id=atom1.id,
                a2_id=atom2.id,
                style=BondStyle.WEDGE,
                stereo=BondStereo.UP,
                stroke_px=stroke,
            )
            item = BondItem(bond, atom1, atom2)
            radians = math.radians(180.0 - angle_deg)
            item.set_wedge_join_neighbors(
                [], [(math.cos(radians), math.sin(radians), stroke, 0.0, 0.0)]
            )
            item.update_positions(atom1, atom2)

            polygons = item.path().toFillPolygons()
            assert len(polygons) == 1
            polygon = polygons[0]
            points = [(polygon[i].x(), polygon[i].y()) for i in range(polygon.size())]
            assert points[0] == pytest.approx(points[-1])
            points = points[:-1]
            assert len(points) >= 3
            assert all(math.isfinite(value) for point in points for value in point)
            area = sum(
                points[i][0] * points[(i + 1) % len(points)][1]
                - points[(i + 1) % len(points)][0] * points[i][1]
                for i in range(len(points))
            )
            assert abs(area) > stroke * stroke
            edges = [
                (points[i], points[(i + 1) % len(points)])
                for i in range(len(points))
            ]
            for i, (a, b) in enumerate(edges):
                for j, (c, d) in enumerate(edges[i + 1 :], start=i + 1):
                    if j == i + 1 or (i == 0 and j == len(edges) - 1):
                        continue
                    assert not _proper_segments_intersect(a, b, c, d)
            assert item.path().boundingRect().width() < length * 1.5
            assert item.path().boundingRect().height() < length * 1.5
    app.processEvents()


def _image_digest(image: QImage) -> str:
    image = image.convertToFormat(QImage.Format.Format_ARGB32)
    ptr = image.bits()
    ptr.setsize(image.sizeInBytes())
    return hashlib.sha256(bytes(ptr)).hexdigest()


@pytest.mark.parametrize(
    "theme,background",
    (("light", "#F1F5F9"), ("dark", "#1B263A")),
)
def test_reference_gallery_is_repeatable_and_theme_surface_is_distinct(theme, background, tmp_path):
    app = QApplication.instance() or QApplication([])
    first = render_gallery(theme, tmp_path / f"bond-gallery-{theme}-a.png")
    second = render_gallery(theme, tmp_path / f"bond-gallery-{theme}-b.png")

    assert first.size() == second.size()
    assert first.pixelColor(0, 0).name() == background.lower()
    assert first.pixelColor(80, 80).name() == "#ffffff"
    assert _image_digest(first) == _image_digest(second)
    app.processEvents()


def test_reference_gallery_can_include_real_atom_numbering(tmp_path):
    app = QApplication.instance() or QApplication([])
    plain = render_gallery("light", tmp_path / "plain.png")
    numbered = render_gallery("light", tmp_path / "numbered.png", numbering=True)

    assert plain.size() == numbered.size()
    assert _image_digest(plain) != _image_digest(numbered)
    app.processEvents()


def test_comparison_sheet_contains_equal_sized_before_after_theme_panels(tmp_path):
    before = tmp_path / "before"
    after = tmp_path / "after"
    before.mkdir()
    after.mkdir()
    for theme in ("light", "dark"):
        for directory, color in ((before, "#FF0000"), (after, "#0000FF")):
            panel = QImage(1680, 1160, QImage.Format.Format_ARGB32_Premultiplied)
            panel.fill(QColor(color))
            assert panel.save(str(directory / f"bond-junctions-{theme}.png"))

    output = tmp_path / "comparison.png"
    compose_comparison(before, after, output)
    sheet = QImage(str(output))

    assert sheet.width() == 3360
    assert sheet.height() == 2 * (1160 + 44)
    assert sheet.pixelColor(10, 44).name() == "#ff0000"
    assert sheet.pixelColor(1690, 44).name() == "#0000ff"
