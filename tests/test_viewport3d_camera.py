from __future__ import annotations

import math

from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
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
