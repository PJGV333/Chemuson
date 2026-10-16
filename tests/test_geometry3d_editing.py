from chemuson.geometry3d.model import CoordinateSet3D
from chemuson.geometry3d.editing import translate_atoms


def test_translate_atoms_changes_only_selected_3d_coordinates():
    original = CoordinateSet3D(
        positions={1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
        source="test",
        method="manual",
    )

    moved = translate_atoms(original, {1, 3}, (0.5, -1.0, 2.0))

    assert moved.positions == {
        1: (0.5, -1.0, 2.0),
        2: (1.0, 0.0, 0.0),
        3: (2.5, -1.0, 2.0),
    }
    assert original.positions[1] == (0.0, 0.0, 0.0)
    assert moved.source == "test"
    assert moved.method == "manual"
