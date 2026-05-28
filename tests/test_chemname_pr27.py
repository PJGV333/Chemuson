"""Pruebas unitarias para la ampliación del motor IUPAC (PR27)."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemname import iupac_name
from chemuson.chemio import rdkit_io
from chemuson.chemio.rdkit_io import smiles_to_molgraph

RDKit_AVAILABLE = rdkit_io._rdkit_available()


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    """Crea una cadena lineal de carbonos y devuelve sus IDs."""
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(length)]
    for i in range(length - 1):
        graph.add_bond(atoms[i].id, atoms[i + 1].id, order=1)
    return [atom.id for atom in atoms]


def build_ring(graph: MolGraph, elements: list[str], double_bonds: set[int]) -> list[int]:
    """Construye un anillo tipo Kekulé con dobles enlaces indicados."""
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    size = len(atoms)
    for i in range(size):
        order = 2 if i in double_bonds else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % size].id, order=order)
    return [atom.id for atom in atoms]


def build_fused_aromatic(graph: MolGraph, second_ring_elements: dict[int, str]) -> list[int]:
    """Construye un sistema fusionado benceno + heterociclo de 5 miembros."""
    elements = ["C"] * 9
    for idx, elem in second_ring_elements.items():
        elements[idx] = elem
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8]
    for ring in (ring1, ring2):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def build_quinoline(graph: MolGraph, n_index: int) -> list[int]:
    """Construye un esqueleto de quinolina/isoquinolina."""
    elements = ["C"] * 10
    elements[n_index] = "N"
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]
    for ring in (ring1, ring2):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def build_naphthalene(graph: MolGraph) -> list[int]:
    """Construye un sistema naftalénico aromático."""
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(10)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8, 9]
    for ring in (ring1, ring2):
        for i in range(len(ring)):
            a1 = atoms[ring[i]].id
            a2 = atoms[ring[(i + 1) % len(ring)]].id
            if graph.find_bond_between(a1, a2) is None:
                graph.add_bond(a1, a2, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def build_benzene_kekule(graph: MolGraph) -> list[int]:
    """Construye benceno en forma Kekulé."""
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(6)]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        graph.add_bond(atoms[i].id, atoms[(i + 1) % 6].id, order=order)
    return [atom.id for atom in atoms]


class ChemNamePR27Test(unittest.TestCase):
    """Casos de ampliación para IUPAC-lite."""

    def test_heterocycles_pyridine_vs_pyridazine(self):
        graph = MolGraph()
        build_ring(graph, ["N", "C", "C", "C", "C", "C"], {0, 2, 4})
        self.assertEqual(iupac_name(graph), "pyridine")

        graph = MolGraph()
        build_ring(graph, ["N", "N", "C", "C", "C", "C"], {0, 2, 4})
        self.assertEqual(iupac_name(graph), "pyridazine")

    def test_fused_quinoline_isoquinoline(self):
        graph = MolGraph()
        build_quinoline(graph, n_index=6)
        self.assertEqual(iupac_name(graph), "quinoline")

        graph = MolGraph()
        build_quinoline(graph, n_index=7)
        self.assertEqual(iupac_name(graph), "isoquinoline")

    def test_fused_benzoxazole_benzimidazole(self):
        graph = MolGraph()
        build_fused_aromatic(graph, {6: "O", 8: "N"})
        self.assertEqual(iupac_name(graph), "benzoxazole")

        graph = MolGraph()
        build_fused_aromatic(graph, {6: "N", 8: "N"})
        self.assertEqual(iupac_name(graph), "benzimidazole")

    def test_spiro_2_5_octane(self):
        graph = MolGraph()
        center = graph.add_atom("C", 0.0, 0.0)
        a1 = graph.add_atom("C", 1.0, 1.0)
        a2 = graph.add_atom("C", 1.0, -1.0)
        b1 = graph.add_atom("C", -1.0, 1.0)
        b2 = graph.add_atom("C", -2.0, 1.0)
        b3 = graph.add_atom("C", -3.0, 0.0)
        b4 = graph.add_atom("C", -2.0, -1.0)
        b5 = graph.add_atom("C", -1.0, -1.0)
        for a, b in (
            (center, a1),
            (a1, a2),
            (a2, center),
            (center, b1),
            (b1, b2),
            (b2, b3),
            (b3, b4),
            (b4, b5),
            (b5, center),
        ):
            graph.add_bond(a.id, b.id, order=1)
        self.assertEqual(iupac_name(graph), "spiro[2.5]octane")

    def test_bicyclo_1_1_1_pentane(self):
        graph = MolGraph()
        a = graph.add_atom("C", 0.0, 0.0)
        b = graph.add_atom("C", 2.0, 0.0)
        c1 = graph.add_atom("C", 1.0, 1.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        c3 = graph.add_atom("C", 1.0, -1.0)
        for x, y in ((a, c1), (c1, b), (a, c2), (c2, b), (a, c3), (c3, b)):
            graph.add_bond(x.id, y.id, order=1)
        self.assertEqual(iupac_name(graph), "bicyclo[1.1.1]pentane")

    def test_bicyclo_4_4_0_decane(self):
        graph = MolGraph()
        atoms = [graph.add_atom("C", float(i), 0.0) for i in range(10)]
        ring1 = [0, 1, 2, 3, 4, 5]
        ring2 = [3, 4, 6, 7, 8, 9]
        for ring in (ring1, ring2):
            for i in range(len(ring)):
                a1 = atoms[ring[i]].id
                a2 = atoms[ring[(i + 1) % len(ring)]].id
                if graph.find_bond_between(a1, a2) is None:
                    graph.add_bond(a1, a2, order=1)
        self.assertEqual(iupac_name(graph), "bicyclo[4.4.0]decane")

    def test_multiple_functional_groups_oxo_hydroxy(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 6)
        hydroxy_o = graph.add_atom("O", 1.0, 1.0, explicit_h=1)
        keto_o = graph.add_atom("O", 2.0, 1.0)
        graph.add_bond(chain[1], hydroxy_o.id, order=1)
        graph.add_bond(chain[2], keto_o.id, order=2)
        self.assertEqual(iupac_name(graph), "2-hydroxyhexan-3-one")

    def test_multiple_functional_groups_amino_oxo_butanoate(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 4)
        carbonyl_o = graph.add_atom("O", 0.0, 1.0)
        ester_o = graph.add_atom("O", 0.0, -1.0)
        ester_methyl = graph.add_atom("C", -1.0, -1.0)
        amino_n = graph.add_atom("N", 1.0, 1.0, explicit_h=1)
        keto_o = graph.add_atom("O", 2.0, 1.0)
        graph.add_bond(chain[0], carbonyl_o.id, order=2)
        graph.add_bond(chain[0], ester_o.id, order=1)
        graph.add_bond(ester_o.id, ester_methyl.id, order=1)
        graph.add_bond(chain[1], amino_n.id, order=1)
        graph.add_bond(chain[2], keto_o.id, order=2)
        self.assertEqual(iupac_name(graph), "2-amino-3-oxobutanoate")

    def test_multiple_unsaturations_linear(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 6)
        for idx in (0, 2, 4):
            graph.find_bond_between(chain[idx], chain[idx + 1]).order = 2
        self.assertEqual(iupac_name(graph), "hexa-1,3,5-triene")

        graph = MolGraph()
        chain = build_linear_chain(graph, 8)
        for idx in (0, 2, 4, 6):
            graph.find_bond_between(chain[idx], chain[idx + 1]).order = 2
        self.assertEqual(iupac_name(graph), "octa-1,3,5,7-tetraene")

    def test_multiple_unsaturations_cyclo(self):
        graph = MolGraph()
        atoms = [graph.add_atom("C", float(i), 0.0) for i in range(8)]
        for i in range(8):
            order = 2 if i in {0, 2, 4, 6} else 1
            graph.add_bond(atoms[i].id, atoms[(i + 1) % 8].id, order=order)
        self.assertEqual(iupac_name(graph), "cyclooct-1,3,5,7-tetraene")

    def test_extended_branched_alkyl_names(self):
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        branch = graph.add_atom("C", -2.0, 0.0)
        m1 = graph.add_atom("C", -3.0, 1.0)
        m2 = graph.add_atom("C", -3.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, branch.id, order=1)
        graph.add_bond(branch.id, m1.id, order=1)
        graph.add_bond(branch.id, m2.id, order=1)
        self.assertEqual(iupac_name(graph), "isobutylbenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        c2 = graph.add_atom("C", -2.0, 0.0)
        branch = graph.add_atom("C", -3.0, 0.0)
        m1 = graph.add_atom("C", -4.0, 1.0)
        m2 = graph.add_atom("C", -4.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, c2.id, order=1)
        graph.add_bond(c2.id, branch.id, order=1)
        graph.add_bond(branch.id, m1.id, order=1)
        graph.add_bond(branch.id, m2.id, order=1)
        self.assertEqual(iupac_name(graph), "isoamylbenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        leaf = graph.add_atom("C", -2.0, 1.0)
        c2 = graph.add_atom("C", -2.0, -1.0)
        c3 = graph.add_atom("C", -3.0, -1.0)
        c4 = graph.add_atom("C", -4.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, leaf.id, order=1)
        graph.add_bond(root.id, c2.id, order=1)
        graph.add_bond(c2.id, c3.id, order=1)
        graph.add_bond(c3.id, c4.id, order=1)
        self.assertEqual(iupac_name(graph), "sec-pentylbenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        l1 = graph.add_atom("C", -2.0, 1.0)
        l2 = graph.add_atom("C", -2.0, -1.0)
        c2 = graph.add_atom("C", -2.0, 0.0)
        c3 = graph.add_atom("C", -3.0, 0.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, l1.id, order=1)
        graph.add_bond(root.id, l2.id, order=1)
        graph.add_bond(root.id, c2.id, order=1)
        graph.add_bond(c2.id, c3.id, order=1)
        self.assertEqual(iupac_name(graph), "tert-pentylbenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        root = graph.add_atom("C", -1.0, 0.0)
        center = graph.add_atom("C", -2.0, 0.0)
        m1 = graph.add_atom("C", -3.0, 1.0)
        m2 = graph.add_atom("C", -3.0, 0.0)
        m3 = graph.add_atom("C", -3.0, -1.0)
        graph.add_bond(ring[0], root.id, order=1)
        graph.add_bond(root.id, center.id, order=1)
        graph.add_bond(center.id, m1.id, order=1)
        graph.add_bond(center.id, m2.id, order=1)
        graph.add_bond(center.id, m3.id, order=1)
        self.assertEqual(iupac_name(graph), "neopentylbenzene")

    def test_extended_alkoxy_benzene(self):
        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        oxygen = graph.add_atom("O", 0.0, 1.0)
        chain = [graph.add_atom("C", float(i), 1.0) for i in range(4)]
        graph.add_bond(ring[0], oxygen.id, order=1)
        graph.add_bond(oxygen.id, chain[0].id, order=1)
        for i in range(3):
            graph.add_bond(chain[i].id, chain[i + 1].id, order=1)
        self.assertEqual(iupac_name(graph), "butoxybenzene")

        graph = MolGraph()
        ring = build_benzene_kekule(graph)
        oxygen = graph.add_atom("O", 0.0, 1.0)
        chain = [graph.add_atom("C", float(i), 1.0) for i in range(5)]
        graph.add_bond(ring[0], oxygen.id, order=1)
        graph.add_bond(oxygen.id, chain[0].id, order=1)
        for i in range(4):
            graph.add_bond(chain[i].id, chain[i + 1].id, order=1)
        self.assertEqual(iupac_name(graph), "pentoxybenzene")

    def test_multiple_substitutions_fused_ring(self):
        graph = MolGraph()
        atoms = build_naphthalene(graph)
        for idx in (5, 0, 1):
            halogen = graph.add_atom("Cl", -1.0, 0.0)
            graph.add_bond(atoms[idx], halogen.id, order=1)
        self.assertEqual(iupac_name(graph), "1,2,3-trichloronaphthalene")

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_stereochemistry_e_2_butene(self):
        graph = smiles_to_molgraph("C/C=C/C")
        self.assertEqual(iupac_name(graph), "(2E)-but-2-ene")

    @unittest.skipIf(not RDKit_AVAILABLE, "RDKit no disponible")
    def test_stereochemistry_lactic_acid(self):
        graph = smiles_to_molgraph("C[C@H](O)C(=O)O")
        self.assertIn(
            iupac_name(graph),
            {
                "(2R)-2-hydroxypropanoic acid",
                "(2S)-2-hydroxypropanoic acid",
            },
        )


if __name__ == "__main__":
    unittest.main()
