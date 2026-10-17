from __future__ import annotations

from chemuson.chemio.rdkit_safe import molgraph_to_inchi_isolated, molgraph_to_smiles_isolated
from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import CoordinateSet3D, OptimizationSettings
from chemuson.geometry3d.model import scene_molecule_from_graph
from chemuson.geometry3d.service import conformer_3d_for_graph


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
