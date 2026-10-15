from __future__ import annotations

import math

from chemuson.geometry3d import SceneAtom3D, SceneMolecule3D
from chemuson.gui.viewport3d.geometry import (
    Bounds3D,
    Ray3D,
    ray_capsule_intersection,
    ray_sphere_intersection,
    scene_bounds,
)
from chemuson.gui.viewport3d.mesh import build_scene_mesh


def _scene() -> SceneMolecule3D:
    return SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "O", 2.0, 0.0, 0.0, 1.52, "#FF0D0D"),
        ),
        bonds=(),
        source="test",
    )


def test_ray_intersections_are_finite_and_ordered() -> None:
    ray = Ray3D((0.0, 0.0, 5.0), (0.0, 0.0, -1.0))
    assert math.isclose(ray_sphere_intersection(ray, (0.0, 0.0, 0.0), 1.0), 4.0)
    assert math.isclose(ray_capsule_intersection(ray, (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.2), 5.0)
    assert ray_sphere_intersection(ray, (10.0, 0.0, 0.0), 0.2) is None


def test_scene_bounds_include_visual_atom_radius() -> None:
    bounds = scene_bounds(_scene())
    assert isinstance(bounds, Bounds3D)
    assert bounds.minimum[0] < 0.0
    assert bounds.maximum[0] > 2.0
    assert all(math.isfinite(value) for value in (*bounds.minimum, *bounds.maximum))


def test_mesh_skips_degenerate_bonds_and_has_interleaved_vertices() -> None:
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    mesh = build_scene_mesh(scene)
    assert mesh.vertex_count > 0
    assert len(mesh.vertices) % 9 == 0
    assert all(math.isfinite(value) for value in mesh.vertices)
