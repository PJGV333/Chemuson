"""Pruebas unitarias PR28: cargas/isótopos/radicales + grupos + heterociclos."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile
from chemuson.chemname import iupac_name
from chemuson.chemname.molview import MolView
from chemuson.chemname.options import NameOptions


def build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(length)]
    for i in range(length - 1):
        graph.add_bond(atoms[i].id, atoms[i + 1].id, order=1)
    return [atom.id for atom in atoms]


def build_aromatic_ring(graph: MolGraph, elements: list[str]) -> list[int]:
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    n = len(atoms)
    for i in range(n):
        graph.add_bond(atoms[i].id, atoms[(i + 1) % n].id, order=1, is_aromatic=True)
    return [atom.id for atom in atoms]


def build_naphthalene(graph: MolGraph) -> list[int]:
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


class ChemNamePR28Test(unittest.TestCase):
    """Cobertura incremental de capacidades PR28 (fases 1-3)."""

    def test_methylammonium_is_not_nd(self):
        graph = MolGraph()
        carbon = graph.add_atom("C", 0.0, 0.0)
        nitrogen = graph.add_atom("N", 1.0, 0.0, charge=1, explicit_h=3)
        graph.add_bond(carbon.id, nitrogen.id, order=1)
        name = iupac_name(graph)
        self.assertNotEqual(name, "N/D")
        self.assertIn("amine", name)

    def test_molview_charge_isotope_radical_accessors(self):
        graph = MolGraph()
        atom = graph.add_atom("O", 0.0, 0.0, charge=-1, isotope=18, radical_electrons=1)
        view = MolView(graph)
        self.assertEqual(view.formal_charge(atom.id), -1)
        self.assertEqual(view.isotope(atom.id), 18)
        self.assertEqual(view.radical_electrons(atom.id), 1)
        self.assertTrue(view.has_radical(atom.id))

    def test_molfile_roundtrip_preserves_isotope_and_radical_without_rdkit_dependency(self):
        graph = MolGraph()
        atom = graph.add_atom("Xx", 0.0, 0.0, charge=-1, isotope=13, radical_electrons=1, is_explicit=True)
        molfile = molgraph_to_molfile(graph)
        restored = molfile_to_molgraph(molfile)
        restored_atom = next(iter(restored.atoms.values()))
        self.assertEqual(restored_atom.formal_charge, -1)
        self.assertEqual(restored_atom.isotope, 13)
        self.assertEqual(getattr(restored_atom, "radical_electrons", 0), 1)

    def test_carboxylate_ethanoate(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        o_dbl = graph.add_atom("O", 1.0, 1.0)
        o_minus = graph.add_atom("O", 1.0, -1.0, charge=-1)
        graph.add_bond(chain[1], o_dbl.id, order=2)
        graph.add_bond(chain[1], o_minus.id, order=1)
        self.assertEqual(iupac_name(graph), "ethanoate")

    def test_dideuterioethane(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        d1 = graph.add_atom("H", -1.0, 1.0, isotope=2, is_explicit=True)
        d2 = graph.add_atom("H", -1.0, -1.0, isotope=2, is_explicit=True)
        graph.add_bond(chain[0], d1.id, order=1)
        graph.add_bond(chain[0], d2.id, order=1)
        self.assertEqual(iupac_name(graph), "1,1-dideuterioethane")

    def test_methoxy_radical_oxyl(self):
        graph = MolGraph()
        carbon = graph.add_atom("C", 0.0, 0.0)
        oxygen = graph.add_atom("O", 1.0, 0.0, radical_electrons=1)
        graph.add_bond(carbon.id, oxygen.id, order=1)
        name = iupac_name(graph)
        self.assertIn("oxyl", name)

    def test_ethanethiol(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        sulfur = graph.add_atom("S", 1.0, 1.0, explicit_h=1)
        graph.add_bond(chain[0], sulfur.id, order=1)
        name = iupac_name(graph)
        self.assertIn("thiol", name)

    def test_dimethyl_sulfoxide_systematic(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        s = graph.add_atom("S", 1.2, 0.0)
        o = graph.add_atom("O", 2.0, 0.8)
        c2 = graph.add_atom("C", 2.4, -0.2)
        graph.add_bond(c1.id, s.id, order=1)
        graph.add_bond(s.id, o.id, order=2)
        graph.add_bond(s.id, c2.id, order=1)
        self.assertEqual(iupac_name(graph), "1-methylsulfinylmethane")

    def test_benzenesulfonic_acid(self):
        graph = MolGraph()
        ring = build_aromatic_ring(graph, ["C", "C", "C", "C", "C", "C"])
        s = graph.add_atom("S", 1.0, 1.0)
        o1 = graph.add_atom("O", 1.8, 1.4)
        o2 = graph.add_atom("O", 1.2, 1.9)
        o3 = graph.add_atom("O", 0.5, 1.4, explicit_h=1)
        graph.add_bond(ring[0], s.id, order=1)
        graph.add_bond(s.id, o1.id, order=2)
        graph.add_bond(s.id, o2.id, order=2)
        graph.add_bond(s.id, o3.id, order=1)
        self.assertEqual(iupac_name(graph), "benzenesulfonic acid")

    def test_azidoethane(self):
        graph = MolGraph()
        chain = build_linear_chain(graph, 2)
        n1 = graph.add_atom("N", 0.0, 1.0)
        n2 = graph.add_atom("N", -1.0, 1.0)
        n3 = graph.add_atom("N", -2.0, 1.0)
        graph.add_bond(chain[0], n1.id, order=1)
        graph.add_bond(n1.id, n2.id, order=2)
        graph.add_bond(n2.id, n3.id, order=1)
        self.assertEqual(iupac_name(graph), "1-azidoethane")

    def test_triazole_and_tetrazole(self):
        graph = MolGraph()
        build_aromatic_ring(graph, ["N", "N", "N", "C", "C"])
        self.assertEqual(iupac_name(graph), "1,2,3-triazole")

        graph = MolGraph()
        build_aromatic_ring(graph, ["N", "N", "N", "N", "C"])
        self.assertEqual(iupac_name(graph), "tetrazole")

    def test_benzotriazole(self):
        graph = MolGraph()
        elements = ["C"] * 9
        elements[6] = "N"
        elements[7] = "N"
        elements[8] = "N"
        atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
        ring1 = [0, 1, 2, 3, 4, 5]
        ring2 = [3, 4, 6, 7, 8]
        for ring in (ring1, ring2):
            for i in range(len(ring)):
                a1 = atoms[ring[i]].id
                a2 = atoms[ring[(i + 1) % len(ring)]].id
                if graph.find_bond_between(a1, a2) is None:
                    graph.add_bond(a1, a2, order=1, is_aromatic=True)
        self.assertEqual(iupac_name(graph), "benzotriazole")

    def test_benzoquinone_and_naphthoquinones(self):
        graph = MolGraph()
        ring = build_aromatic_ring(graph, ["C", "C", "C", "C", "C", "C"])
        o1 = graph.add_atom("O", 0.0, 1.2)
        o2 = graph.add_atom("O", 3.0, 1.2)
        graph.add_bond(ring[0], o1.id, order=2)
        graph.add_bond(ring[3], o2.id, order=2)
        self.assertEqual(iupac_name(graph), "benzene-1,4-dione")

        graph = MolGraph()
        atoms = build_naphthalene(graph)
        o1 = graph.add_atom("O", 6.0, 1.2)
        o2 = graph.add_atom("O", 9.0, 1.2)
        graph.add_bond(atoms[6], o1.id, order=2)
        graph.add_bond(atoms[9], o2.id, order=2)
        self.assertEqual(iupac_name(graph), "naphthalene-1,4-dione")

        graph = MolGraph()
        atoms = build_naphthalene(graph)
        o1 = graph.add_atom("O", 6.0, 1.2)
        o2 = graph.add_atom("O", 7.0, 1.2)
        graph.add_bond(atoms[8], o1.id, order=2)
        graph.add_bond(atoms[9], o2.id, order=2)
        self.assertEqual(iupac_name(graph), "naphthalene-1,2-dione")

    def test_exotic_hetero_flag_phosphabenzene(self):
        graph = MolGraph()
        build_aromatic_ring(graph, ["P", "C", "C", "C", "C", "C"])
        self.assertEqual(iupac_name(graph, NameOptions(enable_exotic_hetero=True)), "phosphabenzene")
        self.assertEqual(iupac_name(graph, NameOptions(enable_exotic_hetero=False)), "N/D")


if __name__ == "__main__":
    unittest.main()
