from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D
from chemuson.geometry3d.editing import rotate_atoms_about_axis, set_dihedral, translate_atoms


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


def _coordinates_with_derived_hydrogens(positions, hydrogen_positions, parents):
    generated_atoms = [
        {"id": hydrogen_id, "element": "H", "coords": list(position)}
        for hydrogen_id, position in hydrogen_positions.items()
    ]
    generated_bonds = [
        {"id": f"bond-{hydrogen_id}", "a1_id": parent_id, "a2_id": hydrogen_id, "order": 1}
        for hydrogen_id, parent_id in parents.items()
    ]
    return CoordinateSet3D(
        positions=positions,
        source="test",
        method="manual",
        metadata={"generated_atoms": generated_atoms, "generated_bonds": generated_bonds},
    )


def _generated_positions(coordinates):
    return {
        item["id"]: tuple(item["coords"])
        for item in coordinates.metadata["generated_atoms"]
    }


def test_translate_atoms_moves_derived_hydrogen_with_parent():
    original = _coordinates_with_derived_hydrogens(
        {1: (0.0, 0.0, 0.0)},
        {"H1": (0.0, 1.0, 0.0)},
        {"H1": 1},
    )

    moved = translate_atoms(original, {1}, (2.0, -1.0, 0.5))

    assert _generated_positions(moved)["H1"] == (2.0, 0.0, 0.5)
    assert _generated_positions(original)["H1"] == (0.0, 1.0, 0.0)


def test_rotate_atoms_about_axis_rotates_derived_hydrogen_with_parent():
    original = _coordinates_with_derived_hydrogens(
        {1: (1.0, 0.0, 0.0)},
        {"H1": (1.0, 1.0, 0.0)},
        {"H1": 1},
    )

    rotated = rotate_atoms_about_axis(original, {1}, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 90.0)

    assert all(
        abs(actual - expected) < 1.0e-9
        for actual, expected in zip(_generated_positions(rotated)["H1"], (-1.0, 1.0, 0.0))
    )
    assert _generated_positions(original)["H1"] == (1.0, 1.0, 0.0)


def test_set_dihedral_rotates_derived_hydrogens_with_moving_fragment():
    graph = MolGraph()
    for atom_id, point in enumerate(
        ((0.0, 1.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 1.0)),
        start=1,
    ):
        graph.add_atom("C", point[0], point[1], atom_id=atom_id)
    for first_id in range(1, 4):
        graph.add_bond(first_id, first_id + 1)
    original = _coordinates_with_derived_hydrogens(
        {
            1: (0.0, 1.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (2.0, 0.0, 0.0),
            4: (3.0, 0.0, 1.0),
        },
        {
            "H1": (0.0, 2.0, 0.0),
            "H2": (1.0, -1.0, 0.0),
            "H3": (2.0, 0.0, 1.0),
            "H4": (3.0, 1.0, 1.0),
        },
        {"H1": 1, "H2": 2, "H3": 3, "H4": 4},
    )

    rotated = set_dihedral(graph, original, (1, 2, 3, 4), 0.0)
    generated = _generated_positions(rotated)

    assert generated["H1"] == (0.0, 2.0, 0.0)
    assert generated["H2"] == (1.0, -1.0, 0.0)
    assert generated["H3"] != (2.0, 0.0, 1.0)
    for hydrogen_id, parent_id in {"H3": 3, "H4": 4}.items():
        parent = rotated.positions[parent_id]
        hydrogen = generated[hydrogen_id]
        assert abs(sum((hydrogen[index] - parent[index]) ** 2 for index in range(3)) ** 0.5 - 1.0) < 1.0e-9
