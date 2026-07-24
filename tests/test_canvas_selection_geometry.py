"""Regresiones directas de la geometría de selección."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
from PyQt6.QtCore import QPointF


GEOMETRY_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "chemuson"
    / "gui"
    / "canvas"
    / "selection_geometry.py"
)
SPEC = importlib.util.spec_from_file_location(
    "chemuson_canvas_selection_geometry",
    GEOMETRY_PATH,
)
assert SPEC is not None and SPEC.loader is not None
GEOMETRY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GEOMETRY)

normalize_custom_stroke = GEOMETRY.normalize_custom_stroke
normalize_label_scale = GEOMETRY.normalize_label_scale
optional_float_equal = GEOMETRY.optional_float_equal
point_equal = GEOMETRY.point_equal
rotate_scene_point = GEOMETRY.rotate_scene_point
scale_point_from_anchor = GEOMETRY.scale_point_from_anchor
signed_angle_delta_deg = GEOMETRY.signed_angle_delta_deg


@pytest.mark.parametrize(
    ("start_deg", "end_deg", "expected"),
    [
        (350.0, 10.0, 20.0),
        (10.0, 350.0, -20.0),
        (0.0, 180.0, -180.0),
        (-90.0, 450.0, -180.0),
    ],
)
def test_signed_angle_delta_preserves_minimal_wrapping(
    start_deg: float,
    end_deg: float,
    expected: float,
) -> None:
    assert signed_angle_delta_deg(start_deg, end_deg) == expected


def test_rotate_scene_point_preserves_center_and_orientation() -> None:
    rotated = rotate_scene_point(QPointF(12.0, 5.0), QPointF(10.0, 5.0), 90.0)

    assert rotated.x() == pytest.approx(10.0)
    assert rotated.y() == pytest.approx(7.0)
    assert rotate_scene_point(QPointF(12.0, 5.0), QPointF(10.0, 5.0), 0.0) == QPointF(12.0, 5.0)


def test_scale_point_from_anchor_preserves_affine_rule() -> None:
    scaled = scale_point_from_anchor(QPointF(2.0, 3.0), QPointF(6.0, 11.0), 0.5)

    assert scaled == QPointF(4.0, 7.0)


def test_optional_float_equality_preserves_none_and_tolerance_rules() -> None:
    assert optional_float_equal(None, None)
    assert not optional_float_equal(None, 1.0)
    assert optional_float_equal(1.0, 1.049, tol=0.05)
    assert not optional_float_equal(1.0, 1.051, tol=0.05)


def test_point_equality_uses_componentwise_tolerance() -> None:
    origin = QPointF(1.0, 2.0)

    assert point_equal(origin, QPointF(1.00009, 1.99991))
    assert not point_equal(origin, QPointF(1.00011, 2.0))


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (0.1, 0.2),
        (0.99, None),
        (1.019, None),
        (1.021, 1.021),
    ],
)
def test_label_scale_normalization_preserves_limits(
    value: float,
    expected: float | None,
) -> None:
    assert normalize_label_scale(value) == expected


@pytest.mark.parametrize(
    ("value", "inherited", "expected"),
    [
        (0.1, 1.0, 0.6),
        (1.03, 1.0, None),
        (1.06, 1.0, 1.06),
        (2.0, 2.0, None),
    ],
)
def test_custom_stroke_normalization_preserves_inheritance_rule(
    value: float,
    inherited: float,
    expected: float | None,
) -> None:
    assert normalize_custom_stroke(value, inherited) == expected
