from __future__ import annotations

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
from chemuson.gui.viewport3d import CameraProjection, CameraState
from chemuson.gui.viewport3d.picking import pick_scene


def test_pick_atom_at_screen_center() -> None:
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    camera = CameraState(target=(0.0, 0.0, 0.0), distance=5.0, yaw=0.0, pitch=0.0)
    hit = pick_scene(scene, camera, 800, 600, 400, 300)
    assert hit is not None
    assert hit.kind == "atom"
    assert hit.primitive_id == 1


def test_pick_bond_when_atom_spheres_do_not_cover_cursor() -> None:
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", -2.0, 0.0, 0.0, 0.1, "#404040"),
            SceneAtom3D(2, "C", 2.0, 0.0, 0.0, 0.1, "#404040"),
        ),
        bonds=(SceneBond3D(7, 1, 2),),
    )
    camera = CameraState(target=(0.0, 0.0, 0.0), distance=5.0, yaw=0.0, pitch=0.0)
    hit = pick_scene(scene, camera, 800, 600, 400, 300)
    assert hit is not None
    assert hit.kind == "bond"
    assert hit.primitive_id == 7


def test_pick_uses_first_surface_depth_when_atom_occludes_bond() -> None:
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(2, "C", -1.0, 0.0, -0.2, 0.1, "#404040"),
            SceneAtom3D(3, "C", 1.0, 0.0, -0.2, 0.1, "#404040"),
        ),
        bonds=(SceneBond3D(9, 2, 3),),
    )
    camera = CameraState(target=(0.0, 0.0, 0.0), distance=5.0, yaw=0.0, pitch=0.0)

    hit = pick_scene(scene, camera, 800, 600, 400, 300)

    assert hit is not None
    assert hit.kind == "atom"
    assert hit.primitive_id == 1


def test_pick_double_bond_outer_lane():
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", -2.0, 0.0, 0.0, 0.01, "#909090"),
            SceneAtom3D(2, "C", 2.0, 0.0, 0.0, 0.01, "#909090"),
            SceneAtom3D(3, "C", -2.0, 1.0, 0.0, 0.01, "#909090"),
        ),
        bonds=(SceneBond3D(7, 1, 2, 2), SceneBond3D(8, 1, 3, 1)),
    )
    camera = CameraState(
        target=(0.0, 0.0, 0.0),
        distance=5.0,
        yaw=0.0,
        pitch=0.0,
        projection=CameraProjection.ORTHOGRAPHIC,
        ortho_scale=10.0,
    )

    hit = pick_scene(scene, camera, 800, 600, 400, 289)

    assert hit is not None
    assert hit.kind == "bond"
    assert hit.primitive_id == 7


def test_picking_radius_tracks_space_filling_representation() -> None:
    scene = SceneMolecule3D(
        atoms=(SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#404040"),),
        bonds=(),
    )
    camera = CameraState(
        target=(0.0, 0.0, 0.0),
        distance=5.0,
        yaw=0.0,
        pitch=0.0,
        projection=CameraProjection.ORTHOGRAPHIC,
        ortho_scale=10.0,
    )

    assert pick_scene(scene, camera, 800, 600, 460, 300) is None
    hit = pick_scene(
        scene,
        camera,
        800,
        600,
        460,
        300,
        representation="space_filling",
    )
    assert hit is not None and hit.primitive_id == 1


def test_pick_preserves_generated_hydrogen_identity_and_parent() -> None:
    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", 8.0, 0.0, 0.0, 1.7, "#404040"),
            SceneAtom3D(
                "genH1",
                "H",
                0.0,
                0.0,
                0.0,
                1.2,
                "#FFFFFF",
                generated=True,
                parent_atom_id=1,
            ),
        ),
        bonds=(SceneBond3D("genB1", 1, "genH1"),),
    )
    camera = CameraState(target=(0.0, 0.0, 0.0), distance=5.0, yaw=0.0, pitch=0.0)

    hit = pick_scene(scene, camera, 800, 600, 400, 300)

    assert hit is not None
    assert hit.kind == "atom"
    assert hit.primitive_id == "genH1"
    assert hit.primitive_id != 1


def test_pick_empty_scene_returns_none() -> None:
    camera = CameraState()
    scene = SceneMolecule3D(atoms=(), bonds=())
    assert pick_scene(scene, camera, 800, 600, 400, 300) is None
