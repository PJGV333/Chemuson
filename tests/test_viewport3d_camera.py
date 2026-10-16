from __future__ import annotations

import math

import pytest

from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
from chemuson.gui.viewport3d.geometry import scene_bounds
from chemuson.gui.viewport3d import CameraProjection, CameraState


def _scene() -> SceneMolecule3D:
    return SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", -2.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "O", 2.0, 0.0, 0.0, 1.52, "#FF0D0D"),
        ),
        bonds=(),
    )


def test_fit_and_navigation_preserve_finite_camera_state() -> None:
    camera = CameraState()
    assert camera.fit_scene(_scene())
    before = camera.target
    camera.orbit(30.0, -20.0, 600.0)
    camera.pan(12.0, -8.0, 600.0)
    camera.zoom(120.0)
    assert camera.target != before
    assert camera.distance > 0.0
    assert all(math.isfinite(value) for value in (*camera.target, camera.distance, camera.yaw, camera.pitch))


def test_orthographic_zoom_changes_scale_not_distance() -> None:
    camera = CameraState(projection=CameraProjection.ORTHOGRAPHIC)
    original_distance = camera.distance
    original_scale = camera.ortho_scale
    camera.zoom(120.0)
    assert camera.distance == original_distance
    assert camera.ortho_scale < original_scale


def test_invalid_viewport_is_rejected_by_ray_parameters() -> None:
    camera = CameraState()
    try:
        camera.ray_parameters(0.0, 500.0)
    except ValueError as exc:
        assert "viewport" in str(exc).lower()
    else:
        raise AssertionError("Se esperaba rechazo de viewport degenerado")


def test_perspective_fit_accounts_for_narrow_viewport_aspect_ratio() -> None:
    scene = _scene()
    camera = CameraState()
    assert camera.fit_scene(scene, viewport_width=120, viewport_height=900)
    bounds = scene_bounds(scene)
    assert bounds is not None
    radius = max(0.05, bounds.radius) * 1.35
    vertical_half_angle = math.radians(camera.fov_y_degrees) / 2.0
    horizontal_half_angle = math.atan(math.tan(vertical_half_angle) * (120 / 900))
    limiting_half_angle = min(vertical_half_angle, horizontal_half_angle)

    assert camera.distance >= radius / math.sin(limiting_half_angle)


def test_orthographic_fit_contains_bounds_in_both_screen_dimensions() -> None:
    scene = _scene()
    camera = CameraState(projection=CameraProjection.ORTHOGRAPHIC)
    assert camera.fit_scene(scene, viewport_width=120, viewport_height=900)
    bounds = scene_bounds(scene)
    assert bounds is not None
    radius = max(0.05, bounds.radius) * 1.35

    assert camera.ortho_scale / 2.0 >= radius
    assert camera.ortho_scale * (120 / 900) / 2.0 + 1.0e-12 >= radius


def test_large_scene_fit_uses_clipping_planes_derived_from_scene() -> None:
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", -1_000_000.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "C", 1_000_000.0, 0.0, 0.0, 1.7, "#404040"),
        ),
        bonds=(),
    )
    camera = CameraState()
    assert camera.fit_scene(scene, viewport_width=1200, viewport_height=800)

    near, far = camera.clip_planes()

    assert camera.distance > 100_000.0
    assert 0.0 < near < camera.distance < far


def test_unrepresentable_viewport_ratio_rejects_fit_without_changing_camera() -> None:
    camera = CameraState(target=(1.0, 2.0, 3.0), distance=12.0, yaw=0.4, pitch=0.2)
    before = (camera.target, camera.distance, camera.fit_radius)

    assert not camera.fit_scene(_scene(), viewport_width=1.0e308, viewport_height=1.0e-308)
    assert (camera.target, camera.distance, camera.fit_radius) == before


def test_screen_world_mapping_matches_perspective_and_orthographic_scale() -> None:
    perspective = CameraState(
        target=(0.0, 0.0, 0.0),
        distance=10.0,
        fov_y_degrees=90.0,
    )
    world = perspective.screen_to_world(
        600.0, 500.0, 1000.0, 1000.0, plane_point=(0.0, 0.0, 0.0)
    )
    assert world == pytest.approx((2.0, 0.0, 0.0), abs=1.0e-8)
    delta = perspective.screen_delta_to_world(
        100.0, 0.0, 1000.0, 1000.0, position=(0.0, 0.0, 0.0)
    )
    assert delta == pytest.approx((2.0, 0.0, 0.0), abs=1.0e-8)

    orthographic = CameraState(
        target=(0.0, 0.0, 0.0),
        ortho_scale=20.0,
        projection=CameraProjection.ORTHOGRAPHIC,
    )
    world = orthographic.screen_to_world(
        600.0, 500.0, 1000.0, 1000.0, plane_point=(0.0, 0.0, 0.0)
    )
    assert world == pytest.approx((2.0, 0.0, 0.0), abs=1.0e-8)
