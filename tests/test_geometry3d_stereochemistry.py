from __future__ import annotations

import math

import pytest

from chemuson.chemio.rdkit_safe import (
    optimize_3d_isolated,
    optimize_3d_stream_isolated,
    smiles_to_molgraph_isolated,
)
from chemuson.core.model import BondStereo, BondStyle, MolGraph
from chemuson.geometry3d.model import OptimizationSettings
from chemuson.geometry3d.service import chemical_graph_hash, conformer_3d_for_graph


@pytest.mark.parametrize(
    ("smiles", "expected"),
    (
        ("F[C@](Cl)(Br)I", "S"),
        ("F[C@@](Cl)(Br)I", "R"),
    ),
)
def test_conformer_preserves_explicit_tetrahedral_cip(smiles, expected):
    graph, error = smiles_to_molgraph_isolated(smiles)
    assert error is None
    assert graph is not None
    center = next(atom for atom in graph.atoms.values() if atom.stereo_cip)
    assert center.stereo_cip == expected

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=40, steps_per_update=10, timeout_s=20.0),
    )

    assert result.ok, result.message
    assert result.metadata.get("stereochemistry") == {
        "atoms": {str(center.id): expected},
        "bonds": {},
    }


@pytest.mark.parametrize(
    ("smiles", "expected"),
    (
        ("F/C=C/Cl", "E"),
        ("F/C=C\\Cl", "Z"),
    ),
)
def test_conformer_preserves_explicit_double_bond_stereo(smiles, expected):
    graph, error = smiles_to_molgraph_isolated(smiles)
    assert error is None
    assert graph is not None
    stereo_bond = next(bond for bond in graph.bonds.values() if bond.stereo_ez)
    assert stereo_bond.stereo_ez == expected

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=40, steps_per_update=10, timeout_s=20.0),
    )

    assert result.ok, result.message
    assert result.metadata.get("stereochemistry") == {
        "atoms": {},
        "bonds": {str(stereo_bond.id): expected},
    }


def test_stereochemical_enantiomers_do_not_share_conformer_cache_key():
    left, left_error = smiles_to_molgraph_isolated("F[C@](Cl)(Br)I")
    right, right_error = smiles_to_molgraph_isolated("F[C@@](Cl)(Br)I")
    assert left_error is None
    assert right_error is None
    assert left is not None
    assert right is not None

    assert chemical_graph_hash(left) != chemical_graph_hash(right)


def test_unspecified_stereochemistry_is_not_reported_as_preserved():
    graph, error = smiles_to_molgraph_isolated("FCCCl")
    assert error is None
    assert graph is not None

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=20, steps_per_update=10, timeout_s=20.0),
    )

    assert result.ok, result.message
    assert result.metadata.get("stereochemistry") == {"atoms": {}, "bonds": {}}


@pytest.mark.parametrize(
    ("smiles", "expected"),
    (
        ("F[C@](Cl)(Br)I", "S"),
        ("F/C=C/Cl", "E"),
    ),
)
def test_forcefield_optimization_preserves_stereo_and_document_atom_ids(smiles, expected):
    graph, error = smiles_to_molgraph_isolated(smiles)
    assert error is None
    assert graph is not None

    positions, metadata, error = optimize_3d_isolated(
        graph, max_iters=30, steps_per_update=10, timeout_s=20.0
    )

    assert error is None
    assert positions is not None
    assert set(positions) == set(graph.atoms)
    labels = metadata.get("stereochemistry")
    assert isinstance(labels, dict)
    if expected in {"R", "S"}:
        center = next(atom for atom in graph.atoms.values() if atom.stereo_cip)
        assert labels == {"atoms": {str(center.id): expected}, "bonds": {}}
    else:
        stereo_bond = next(bond for bond in graph.bonds.values() if bond.stereo_ez)
        assert labels == {"atoms": {}, "bonds": {str(stereo_bond.id): expected}}


def test_forcefield_rejects_inverted_initial_tetrahedral_coordinates():
    graph, error = smiles_to_molgraph_isolated("F[C@](Cl)(Br)I")
    assert error is None
    assert graph is not None
    conformer = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=30, steps_per_update=10, timeout_s=20.0),
    )
    assert conformer.ok, conformer.message
    center = next(atom for atom in graph.atoms.values() if atom.stereo_cip)
    neighbors = sorted(
        bond.a2_id if bond.a1_id == center.id else bond.a1_id
        for bond in graph.bonds.values()
        if center.id in (bond.a1_id, bond.a2_id)
    )
    inverted_coordinates = dict(conformer.atom_positions)
    inverted_coordinates[neighbors[0]], inverted_coordinates[neighbors[1]] = (
        inverted_coordinates[neighbors[1]],
        inverted_coordinates[neighbors[0]],
    )

    positions, _metadata, error = optimize_3d_isolated(
        graph,
        coordinates=inverted_coordinates,
        max_iters=10,
        steps_per_update=5,
        timeout_s=20.0,
    )

    assert positions is None
    assert error is not None
    assert error.startswith(f"stereochemistry_inversion:atom:{center.id}:")


def test_forcefield_rejects_inverted_initial_ez_coordinates():
    graph, error = smiles_to_molgraph_isolated("F/C=C/Cl")
    assert error is None
    assert graph is not None
    stereo_bond = next(bond for bond in graph.bonds.values() if bond.stereo_ez)
    conformer = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=30, steps_per_update=10, timeout_s=20.0),
    )
    assert conformer.ok, conformer.message

    positions = dict(conformer.atom_positions)
    axis_start = positions[stereo_bond.a1_id]
    axis_end = positions[stereo_bond.a2_id]
    axis_vector = tuple(axis_end[index] - axis_start[index] for index in range(3))
    axis_length = math.sqrt(sum(value * value for value in axis_vector))
    axis = (
        axis_vector[0] / axis_length,
        axis_vector[1] / axis_length,
        axis_vector[2] / axis_length,
    )
    fluorine_id = next(atom.id for atom in graph.atoms.values() if atom.element == "F")
    relative = tuple(positions[fluorine_id][index] - axis_start[index] for index in range(3))
    projection = sum(relative[index] * axis[index] for index in range(3))
    positions[fluorine_id] = tuple(
        axis_start[index] + 2.0 * projection * axis[index] - relative[index]
        for index in range(3)
    )

    optimized, _metadata, error = optimize_3d_isolated(
        graph,
        coordinates=positions,
        max_iters=10,
        steps_per_update=5,
        timeout_s=20.0,
    )

    assert optimized is None
    assert error is not None
    assert error.startswith(f"stereochemistry_inversion:bond:{stereo_bond.id}:")


def test_live_forcefield_frames_report_only_preserved_stereo():
    graph, error = smiles_to_molgraph_isolated("F/C=C\\Cl")
    assert error is None
    assert graph is not None
    stereo_bond = next(bond for bond in graph.bonds.values() if bond.stereo_ez)
    expected = {"atoms": {}, "bonds": {str(stereo_bond.id): "Z"}}

    events = list(
        optimize_3d_stream_isolated(
            graph, max_iters=20, steps_per_update=10, timeout_s=20.0
        )
    )

    frames = [event for event in events if event.get("event") == "frame"]
    assert frames
    assert all(frame.get("stereochemistry") == expected for frame in frames)
    assert events[-1].get("event") == "done"
    assert events[-1].get("metadata", {}).get("stereochemistry") == expected


def test_manual_wedge_without_cip_metadata_is_interpreted_from_2d_graph():
    graph = MolGraph()
    fluorine = graph.add_atom("F", 1.5, 0.0).id
    center = graph.add_atom("C", 0.0, 0.0).id
    chlorine = graph.add_atom("Cl", -1.5, 0.0).id
    bromine = graph.add_atom("Br", 0.0, 1.5).id
    iodine = graph.add_atom("I", 0.0, -1.5).id
    graph.add_bond(center, fluorine, style=BondStyle.HASHED, stereo=BondStereo.DOWN)
    graph.add_bond(center, chlorine)
    graph.add_bond(center, bromine)
    graph.add_bond(center, iodine)

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=30, steps_per_update=10, timeout_s=20.0),
    )

    assert result.ok, result.message
    assert result.metadata.get("stereochemistry") == {
        "atoms": {str(center): "S"},
        "bonds": {},
    }
