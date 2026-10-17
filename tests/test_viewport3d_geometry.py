from __future__ import annotations

import math

import pytest

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
from chemuson.gui.viewport3d.geometry import (
    ATOM_RADIUS_SCALE,
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
    assert ray_capsule_intersection(ray, (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.2) == pytest.approx(4.8)
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


def test_mesh_preserves_single_double_triple_and_aromatic_bond_rendering() -> None:
    atoms = (
        SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 0.0, "#404040"),
        SceneAtom3D(2, "C", 1.4, 0.0, 0.0, 0.0, "#404040"),
        SceneAtom3D(3, "C", 0.0, 1.0, 0.0, 0.0, "#404040"),
    )

    def mesh_for(order: int, *, aromatic: bool = False):
        scene = SceneMolecule3D(
            atoms=atoms,
            bonds=(SceneBond3D(1, 1, 2, order, aromatic),),
        )
        return build_scene_mesh(scene)

    single = mesh_for(1)
    double = mesh_for(2)
    triple = mesh_for(3)
    aromatic = mesh_for(1, aromatic=True)

    assert single.vertex_count < double.vertex_count < triple.vertex_count
    assert single.vertex_count < aromatic.vertex_count
    assert aromatic.vertices != double.vertices
    assert all(
        math.isfinite(value)
        for value in (*single.vertices, *double.vertices, *triple.vertices, *aromatic.vertices)
    )


def test_ball_and_stick_radii_leave_room_to_see_carbon_hydrogen_bonds() -> None:
    carbon_radius = 1.7 * ATOM_RADIUS_SCALE
    hydrogen_radius = 1.2 * ATOM_RADIUS_SCALE

    assert carbon_radius + hydrogen_radius < 1.09


def test_capsule_pick_returns_first_lateral_surface_contact() -> None:
    ray = Ray3D((0.0, 2.0, 0.0), (0.0, -1.0, 0.0))

    distance = ray_capsule_intersection(ray, (-2.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.25)

    assert distance == pytest.approx(1.75)


def test_capsule_pick_from_inside_exits_at_the_first_surface() -> None:
    ray = Ray3D((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))

    distance = ray_capsule_intersection(ray, (-2.0, 0.0, 0.0), (2.0, 0.0, 0.0), 0.25)

    assert distance == pytest.approx(0.25)


def test_capsule_pick_from_inside_along_axis_exits_through_outer_cap() -> None:
    ray = Ray3D((0.5, 0.0, 0.0), (1.0, 0.0, 0.0))

    distance = ray_capsule_intersection(ray, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.25)

    assert distance == pytest.approx(0.75)


def test_capsule_pick_hits_only_the_outward_hemisphere_of_an_endcap() -> None:
    ray = Ray3D((2.0, 0.0, 0.0), (-1.0, 0.0, 0.0))

    distance = ray_capsule_intersection(ray, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.25)

    assert distance == pytest.approx(0.75)


def test_capsule_pick_handles_tangent_and_nearly_parallel_rays() -> None:
    tangent = Ray3D((-2.0, 0.25, 0.0), (1.0, 0.0, 0.0))
    miss = Ray3D((-2.0, 0.251, 0.0), (1.0, 1.0e-10, 0.0))

    assert ray_capsule_intersection(tangent, (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.25) == pytest.approx(1.0)
    assert ray_capsule_intersection(miss, (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.25) is None


def test_capsule_pick_handles_short_and_degenerate_segments() -> None:
    ray = Ray3D((0.0, 1.0, 0.0), (0.0, -1.0, 0.0))

    short = ray_capsule_intersection(ray, (-0.01, 0.0, 0.0), (0.01, 0.0, 0.0), 0.25)
    degenerate = ray_capsule_intersection(ray, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.25)

    assert short == pytest.approx(0.75)
    assert degenerate == pytest.approx(0.75)
