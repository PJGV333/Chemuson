from __future__ import annotations

import math
from pathlib import Path

from chemuson.chemio.rdkit_io import molfile_to_molgraph
from chemuson.chemio.rdkit_safe import molgraph_to_inchi_isolated, molgraph_to_smiles_isolated
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D, OptimizationSettings
from chemuson.geometry3d.model import scene_molecule_from_graph
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


def test_local_tetrandrina_template_reaches_scene_mesh_without_inventing_stereo() -> None:
    template_path = Path(__file__).parent / "data" / "tetrandrina_local_template.mol"
    graph = molfile_to_molgraph(template_path.read_text(encoding="utf-8"))
    assert len(graph.atoms) == 43
    assert len(graph.bonds) == 48
    assert sum(atom.element == "CH3" for atom in graph.atoms.values()) == 2
    assert sum(atom.element == "OCH3" for atom in graph.atoms.values()) == 4

    original_atoms = {
        atom_id: (atom.element, atom.charge, atom.isotope, atom.stereo_cip)
        for atom_id, atom in graph.atoms.items()
    }
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
    coordinates = result.coordinate_set
    assert coordinates is not None
    diagnostics = coordinates.metadata.get("diagnostics")
    assert isinstance(diagnostics, dict)
    assert diagnostics.get("stage") == "done"
    stereo_diagnostics = diagnostics.get("stereochemistry_diagnostics")
    assert isinstance(stereo_diagnostics, dict)
    assert stereo_diagnostics.get("specified_atoms") == {}
    assert stereo_diagnostics.get("specified_bonds") == {}
    potential_atoms = stereo_diagnostics.get("potential_atom_centers")
    potential_bonds = stereo_diagnostics.get("potential_bond_stereo")
    assert isinstance(potential_atoms, list) and all(
        not item.get("specified", True) for item in potential_atoms if isinstance(item, dict)
    )
    assert isinstance(potential_bonds, list) and all(
        not item.get("specified", True) for item in potential_bonds if isinstance(item, dict)
    )
    minimization = diagnostics.get("minimization")
    assert isinstance(minimization, dict)
    assert minimization.get("converged") is False
    assert coordinates.energy is None
    assert result.method == "unoptimized"
    template_normalization = diagnostics.get("template_normalization")
    assert isinstance(template_normalization, dict)
    assert template_normalization.get("aromatic_bridge_bond_ids") == [9]

    scene = scene_molecule_from_graph(graph, coordinates)
    scene_atoms = {atom.id: atom for atom in scene.atoms}
    assert scene_atoms[12].element == "C"
    assert scene_atoms[20].element == "O"
    assert len(
        [atom for atom in scene.atoms if atom.generated and atom.element == "C"]
    ) == 4
    assert len(scene.atoms) >= len(graph.atoms) + 4
    mesh = build_scene_mesh(scene)
    assert mesh.vertex_count > 0
    assert all(math.isfinite(value) for value in mesh.vertices)
    assert {
        atom_id: (atom.element, atom.charge, atom.isotope, atom.stereo_cip)
        for atom_id, atom in graph.atoms.items()
    } == original_atoms
    assert {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic, bond.stereo_ez)
        for bond_id, bond in graph.bonds.items()
    } == original_bonds
