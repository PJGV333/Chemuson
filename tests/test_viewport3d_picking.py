from __future__ import annotations

from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D
from chemuson.gui.viewport3d import CameraState
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


def test_pick_empty_scene_returns_none() -> None:
    camera = CameraState()
    scene = SceneMolecule3D(atoms=(), bonds=())
    assert pick_scene(scene, camera, 800, 600, 400, 300) is None
