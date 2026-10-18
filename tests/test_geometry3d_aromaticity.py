from __future__ import annotations

import math
from pathlib import Path

import pytest

from chemuson.chemio.rdkit_io import (
    expand_abbreviations_for_calculation,
    molfile_to_molgraph,
)
from chemuson.chemio.rdkit_safe import (
    molecular_descriptors_isolated,
    molgraph_to_inchi_isolated,
    molgraph_to_smiles_isolated,
    smiles_to_molgraph_isolated,
)
from chemuson.core.model import MolGraph
from chemuson.geometry3d import rdkit_backend
from chemuson.geometry3d.model import (
    CoordinateSet3D,
    OptimizationSettings,
    scene_molecule_from_graph,
)
from chemuson.geometry3d.service import conformer_3d_for_graph
from chemuson.gui.viewport3d.mesh import build_scene_mesh


_RING_COORDS = (
    (1.0, 0.0),
    (0.5, 0.866),
    (-0.5, 0.866),
    (-1.0, 0.0),
    (-0.5, -0.866),
    (0.5, -0.866),
)


def _benzene_graph(*, aromatic: bool) -> MolGraph:
    graph = MolGraph()
    atoms = [graph.add_atom("C", x, y).id for x, y in _RING_COORDS]
    for index, atom_id in enumerate(atoms):
        graph.add_bond(
            atom_id,
            atoms[(index + 1) % len(atoms)],
            order=1 if aromatic or index % 2 == 1 else 2,
            is_aromatic=aromatic,
        )
    return graph


def test_aromatic_and_kekule_benzene_have_the_same_canonical_chemistry():
    aromatic = _benzene_graph(aromatic=True)
    kekule = _benzene_graph(aromatic=False)

    aromatic_smiles, aromatic_smiles_error = molgraph_to_smiles_isolated(aromatic)
    kekule_smiles, kekule_smiles_error = molgraph_to_smiles_isolated(kekule)
    aromatic_inchi, aromatic_inchi_error = molgraph_to_inchi_isolated(aromatic)
    kekule_inchi, kekule_inchi_error = molgraph_to_inchi_isolated(kekule)

    assert aromatic_smiles_error is None
    assert kekule_smiles_error is None
    assert aromatic_smiles == "c1ccccc1"
    assert kekule_smiles == "C1=CC=CC=C1"
    assert aromatic_inchi_error is None
    assert kekule_inchi_error is None
    assert aromatic_inchi == kekule_inchi




def test_scene_molecule_retains_aromatic_flags_and_kekule_orders():
    for graph in (_benzene_graph(aromatic=True), _benzene_graph(aromatic=False)):
        coordinates = CoordinateSet3D(
            {
                atom.id: (atom.x, atom.y, 0.0)
                for atom in graph.atoms.values()
            },
            source="test",
            method="manual",
        )
        scene = scene_molecule_from_graph(graph, coordinates)

        assert [(bond.order, bond.is_aromatic) for bond in scene.bonds] == [
            (bond.order, bond.is_aromatic) for bond in graph.bonds.values()
        ]


def test_aromatic_and_kekule_benzene_both_generate_3d_coordinates():
    for graph in (_benzene_graph(aromatic=True), _benzene_graph(aromatic=False)):
        result = conformer_3d_for_graph(
            graph,
            force=True,
            settings=OptimizationSettings(max_iters=40, steps_per_update=10, timeout_s=20.0),
        )

        assert result.ok, result.message
        assert set(result.atom_positions) == set(graph.atoms)
        assert len(result.atom_positions) == 6


_PUBCHEM_TETRANDRINA_CONNECTIVITY_SMILES = (
    "CN1CCC2=CC(=C3C=C2C1CC4=CC=C(C=C4)OC5="
    "C(C=CC(=C5)CC6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"
)


def _connected_component_count(graph: MolGraph) -> int:
    adjacency = {atom_id: set() for atom_id in graph.atoms}
    for bond in graph.bonds.values():
        adjacency[bond.a1_id].add(bond.a2_id)
        adjacency[bond.a2_id].add(bond.a1_id)
    unseen = set(adjacency)
    components = 0
    while unseen:
        components += 1
        pending = [unseen.pop()]
        while pending:
            current = pending.pop()
            neighbors = adjacency[current] & unseen
            unseen.difference_update(neighbors)
            pending.extend(neighbors)
    return components


def test_acyclic_aromatic_bridge_is_rejected_without_mutating_graph() -> None:
    graph = MolGraph()
    first = graph.add_atom("C", 0.0, 0.0).id
    middle = graph.add_atom("C", 1.0, 0.0).id
    last = graph.add_atom("C", 2.0, 0.0).id
    aromatic_bridge = graph.add_bond(first, middle, order=1, is_aromatic=True)
    graph.add_bond(middle, last, order=1)

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=20, timeout_s=10.0),
    )

    assert not result.ok
    diagnostics = result.metadata.get("diagnostics")
    assert isinstance(diagnostics, dict)
    assert diagnostics.get("error_code") == "invalid_aromatic_bridge"
    assert diagnostics.get("stage") == "calculation_graph_validation"
    assert diagnostics.get("aromatic_bridge_bond_ids") == [aromatic_bridge.id]
    assert graph.bonds[aromatic_bridge.id].is_aromatic is True


def test_optimizer_entry_points_reject_aromatic_bridges_consistently() -> None:
    graph = MolGraph()
    first = graph.add_atom("C", 0.0, 0.0).id
    middle = graph.add_atom("C", 1.0, 0.0).id
    last = graph.add_atom("C", 2.0, 0.0).id
    aromatic_bridge = graph.add_bond(first, middle, order=1, is_aromatic=True)
    graph.add_bond(middle, last, order=1)
    settings = OptimizationSettings(max_iters=20, steps_per_update=10, timeout_s=10.0)

    for result in (
        rdkit_backend.optimize(graph, settings=settings),
        rdkit_backend.optimize_stream(graph, settings=settings),
    ):
        assert not result.ok
        diagnostics = result.metadata.get("diagnostics")
        assert isinstance(diagnostics, dict)
        assert diagnostics.get("error_code") == "invalid_aromatic_bridge"
        assert diagnostics.get("stage") == "calculation_graph_validation"
        assert diagnostics.get("aromatic_bridge_bond_ids") == [aromatic_bridge.id]

    with pytest.raises(RuntimeError, match="invalid_aromatic_bridge") as error:
        list(rdkit_backend.optimize_iter(graph, settings=settings))
    assert getattr(error.value, "diagnostics", {}).get("stage") == "calculation_graph_validation"
    assert graph.bonds[aromatic_bridge.id].is_aromatic is True


def test_historical_tetrandrina_fixture_is_rejected_without_hidden_repair() -> None:
    template_path = Path(__file__).parent / "data" / "tetrandrina_local_template.mol"
    graph = molfile_to_molgraph(template_path.read_text(encoding="utf-8"))
    original_atoms = {
        atom_id: (atom.element, atom.x, atom.y, atom.stereo_cip)
        for atom_id, atom in graph.atoms.items()
    }
    original_bonds = {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic, bond.stereo_ez)
        for bond_id, bond in graph.bonds.items()
    }

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=40, timeout_s=10.0, seed=173),
    )

    assert not result.ok
    diagnostics = result.metadata.get("diagnostics")
    assert isinstance(diagnostics, dict)
    assert diagnostics.get("error_code") == "invalid_aromatic_bridge"
    assert diagnostics.get("stage") == "calculation_graph_validation"
    assert diagnostics.get("aromatic_bridge_bond_ids") == [9]
    assert {atom_id: (atom.element, atom.x, atom.y, atom.stereo_cip) for atom_id, atom in graph.atoms.items()} == original_atoms
    assert {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic, bond.stereo_ez)
        for bond_id, bond in graph.bonds.items()
    } == original_bonds


def test_corrected_tetrandrina_fixture_matches_pubchem_connectivity_and_expands() -> None:
    source_path = Path(__file__).parent / "data" / "tetrandrina_local_template.mol"
    corrected_path = Path(__file__).parent / "data" / "tetrandrina_local_template_repaired.mol"
    source_graph = molfile_to_molgraph(source_path.read_text(encoding="utf-8"))
    graph = molfile_to_molgraph(corrected_path.read_text(encoding="utf-8"))

    assert len(graph.atoms) == 42
    assert len(graph.bonds) == 48
    assert sum(atom.element == "CH3" for atom in graph.atoms.values()) == 2
    assert sum(atom.element == "OCH3" for atom in graph.atoms.values()) == 4
    closure = graph.find_bond_between(5, 6)
    assert closure is not None
    assert closure.order == 1 and closure.is_aromatic

    old_to_new = {
        atom_id: atom_id if atom_id < 11 else atom_id - 1
        for atom_id in source_graph.atoms
        if atom_id != 11
    }
    assert set(old_to_new.values()) == set(graph.atoms)
    for old_id, new_id in old_to_new.items():
        old_atom = source_graph.atoms[old_id]
        new_atom = graph.atoms[new_id]
        assert (
            old_atom.element,
            old_atom.charge,
            old_atom.isotope,
            old_atom.stereo_cip,
            old_atom.explicit_h,
        ) == (
            new_atom.element,
            new_atom.charge,
            new_atom.isotope,
            new_atom.stereo_cip,
            new_atom.explicit_h,
        )
        assert math.isclose(old_atom.x, new_atom.x, abs_tol=1.0e-3)
        assert math.isclose(old_atom.y, new_atom.y, abs_tol=1.0e-3)
    for old_bond_id, old_bond in source_graph.bonds.items():
        if old_bond_id == 9:
            continue
        new_bond = graph.find_bond_between(old_to_new[old_bond.a1_id], old_to_new[old_bond.a2_id])
        assert new_bond is not None
        assert (new_bond.order, new_bond.is_aromatic, new_bond.stereo_ez) == (
            old_bond.order,
            old_bond.is_aromatic,
            old_bond.stereo_ez,
        )

    expanded = expand_abbreviations_for_calculation(graph)
    assert len(expanded.atoms) == 46
    assert len(expanded.bonds) == 52
    element_counts: dict[str, int] = {}
    for atom in expanded.atoms.values():
        element_counts[atom.element] = element_counts.get(atom.element, 0) + 1
    assert element_counts == {"C": 38, "N": 2, "O": 6}
    assert _connected_component_count(expanded) == 1
    assert len(expanded.bonds) - len(expanded.atoms) + 1 == 7
    aromatic_bonds = [bond for bond in expanded.bonds.values() if bond.is_aromatic]
    aromatic_atoms = {
        atom_id
        for bond in aromatic_bonds
        for atom_id in (bond.a1_id, bond.a2_id)
    }
    assert len(aromatic_bonds) == 24
    assert len(aromatic_atoms) == 24
    endpoint_pairs = [
        tuple(sorted((bond.a1_id, bond.a2_id)))
        for bond in expanded.bonds.values()
    ]
    assert len(endpoint_pairs) == len(set(endpoint_pairs))
    assert all(first != second for first, second in endpoint_pairs)
    assert expanded.validate_detailed() == {}
    assert all(atom.stereo_cip is None for atom in graph.atoms.values())
    assert all(bond.stereo_ez is None for bond in graph.bonds.values())

    inchi, inchi_error = molgraph_to_inchi_isolated(expanded)
    reference_graph, reference_error = smiles_to_molgraph_isolated(
        _PUBCHEM_TETRANDRINA_CONNECTIVITY_SMILES
    )
    assert inchi_error is None and inchi is not None
    assert reference_error is None and reference_graph is not None
    reference_inchi, reference_inchi_error = molgraph_to_inchi_isolated(reference_graph)
    assert reference_inchi_error is None and reference_inchi == inchi
    assert inchi.split("/")[1] == "C38H42N2O6"
    assert len(expanded.atoms) + 42 == 88

    descriptors, descriptor_error = molecular_descriptors_isolated(expanded)
    assert descriptor_error is None and descriptors is not None
    assert descriptors.get("heavy_atoms") == 46
    smiles, smiles_error = molgraph_to_smiles_isolated(expanded)
    assert smiles_error is None and smiles
    smiles_graph, smiles_graph_error = smiles_to_molgraph_isolated(smiles)
    assert smiles_graph_error is None and smiles_graph is not None
    smiles_inchi, smiles_inchi_error = molgraph_to_inchi_isolated(smiles_graph)
    assert smiles_inchi_error is None and smiles_inchi == inchi


def test_corrected_tetrandrina_reaches_worker_scene_and_reports_mmff_status() -> None:
    template_path = Path(__file__).parent / "data" / "tetrandrina_local_template_repaired.mol"
    graph = molfile_to_molgraph(template_path.read_text(encoding="utf-8"))
    original_bonds = {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic, bond.stereo_ez)
        for bond_id, bond in graph.bonds.items()
    }

    result = conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(max_iters=80, steps_per_update=10, timeout_s=30.0, seed=173),
    )

    assert result.ok, result.message
    assert len(result.atom_positions) == 46
    assert all(math.isfinite(value) for point in result.atom_positions.values() for value in point)
    diagnostics = result.metadata.get("diagnostics")
    assert isinstance(diagnostics, dict)
    assert diagnostics.get("stage") == "done"
    stereo = diagnostics.get("stereochemistry_diagnostics")
    assert isinstance(stereo, dict)
    assert stereo.get("specified_atoms") == {}
    assert stereo.get("specified_bonds") == {}
    minimization = diagnostics.get("minimization")
    assert isinstance(minimization, dict)
    assert minimization.get("requested_forcefield") == "MMFF94"
    assert minimization.get("selected_forcefield") == "mmff94"
    assert isinstance(minimization.get("converged"), bool)
    if minimization["converged"]:
        assert result.method != "unoptimized"
        assert result.energy is not None and math.isfinite(result.energy)
    else:
        assert minimization.get("detail") == "minimization_not_converged"
        assert result.method == "unoptimized"
        assert result.energy is None

    coordinates = result.coordinate_set
    assert coordinates is not None
    scene = scene_molecule_from_graph(graph, coordinates)
    assert len(scene.atoms) == 88
    assert len(scene.bonds) == 94
    mesh = build_scene_mesh(scene)
    assert mesh.vertex_count > 0
    assert all(math.isfinite(value) for value in mesh.vertices)
    assert {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic, bond.stereo_ez)
        for bond_id, bond in graph.bonds.items()
    } == original_bonds
