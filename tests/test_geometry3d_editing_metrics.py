import pytest

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D
from chemuson.geometry3d.editing import (
    angle_degrees,
    bond_is_bridge,
    dihedral_degrees,
    distance,
    set_bond_length,
    set_bond_angle,
    set_dihedral,
)


def _butane():
    graph = MolGraph()
    atoms = [graph.add_atom("C", index * 40.0, 0.0) for index in range(4)]
    for first, second in zip(atoms, atoms[1:]):
        graph.add_bond(first.id, second.id, order=1)
    coordinates = CoordinateSet3D(
        {
            atoms[0].id: (0.0, 1.0, 0.0),
            atoms[1].id: (0.0, 0.0, 0.0),
            atoms[2].id: (1.0, 0.0, 0.0),
            atoms[3].id: (1.0, 1.0, 1.0),
        },
        source="test",
        method="manual",
    )
    return graph, [atom.id for atom in atoms], coordinates


def test_set_bond_length_moves_the_entire_bridge_side():
    graph, atom_ids, coordinates = _butane()
    bond = next(b for b in graph.bonds.values() if {b.a1_id, b.a2_id} == {atom_ids[1], atom_ids[2]})

    moved = set_bond_length(graph, coordinates, bond.id, 2.0, moving_endpoint_id=atom_ids[2])

    assert distance(moved, atom_ids[1], atom_ids[2]) == pytest.approx(2.0)
    assert distance(moved, atom_ids[2], atom_ids[3]) == pytest.approx(
        distance(coordinates, atom_ids[2], atom_ids[3])
    )
    assert moved.positions[atom_ids[0]] == coordinates.positions[atom_ids[0]]
    assert moved.positions[atom_ids[1]] == coordinates.positions[atom_ids[1]]


def test_set_angle_rotates_branch_and_reaches_target():
    graph, atom_ids, coordinates = _butane()
    result = set_bond_angle(graph, coordinates, (atom_ids[0], atom_ids[1], atom_ids[2]), 120.0)

    assert angle_degrees(result, atom_ids[0], atom_ids[1], atom_ids[2]) == pytest.approx(120.0)
    assert result.positions[atom_ids[0]] == coordinates.positions[atom_ids[0]]
    assert result.positions[atom_ids[1]] == coordinates.positions[atom_ids[1]]


def test_set_dihedral_reaches_target_without_moving_fixed_side():
    graph, atom_ids, coordinates = _butane()
    assert bond_is_bridge(graph, next(b.id for b in graph.bonds.values() if b.a1_id == atom_ids[1] and b.a2_id == atom_ids[2]))

    moved = set_dihedral(graph, coordinates, atom_ids, 60.0)

    assert dihedral_degrees(moved, *atom_ids) == pytest.approx(60.0, abs=1e-6)
    assert moved.positions[atom_ids[0]] == coordinates.positions[atom_ids[0]]
    assert moved.positions[atom_ids[1]] == coordinates.positions[atom_ids[1]]


def test_ring_bond_is_not_a_rotatable_independent_fragment():
    graph = MolGraph()
    atoms = [graph.add_atom("C", index * 40.0, 0.0) for index in range(4)]
    for index in range(4):
        graph.add_bond(atoms[index].id, atoms[(index + 1) % 4].id, order=1)
    coordinates = CoordinateSet3D(
        {
            atoms[0].id: (0.0, 0.0, 0.0),
            atoms[1].id: (1.0, 0.0, 0.0),
            atoms[2].id: (1.0, 1.0, 0.0),
            atoms[3].id: (0.0, 1.0, 0.0),
        },
        source="test",
        method="manual",
    )
    ring_bond = next(b for b in graph.bonds.values() if {b.a1_id, b.a2_id} == {atoms[1].id, atoms[2].id})

    assert not bond_is_bridge(graph, ring_bond.id)
    with pytest.raises(ValueError, match="bond_not_bridge"):
        set_dihedral(graph, coordinates, tuple(atom.id for atom in atoms), 90.0)
