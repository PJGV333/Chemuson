import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from core.model import MolGraph
from chemname import iupac_name
from chemname.options import NameOptions
from chemname.template import load_template


def add_aromatic_ring(graph: MolGraph, atoms: list, ring: list[int]) -> None:
    for i in range(len(ring)):
        a1 = atoms[ring[i]].id
        a2 = atoms[ring[(i + 1) % len(ring)]].id
        if graph.find_bond_between(a1, a2) is None:
            graph.add_bond(a1, a2, order=1, is_aromatic=True)


def build_benzofuran(graph: MolGraph) -> list[int]:
    elements = ["C"] * 9
    elements[6] = "O"
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8]
    add_aromatic_ring(graph, atoms, ring1)
    add_aromatic_ring(graph, atoms, ring2)
    return [atom.id for atom in atoms]


def build_benzothiophene(graph: MolGraph) -> list[int]:
    elements = ["C"] * 9
    elements[6] = "S"
    atoms = [graph.add_atom(elem, float(i), 0.0) for i, elem in enumerate(elements)]
    ring1 = [0, 1, 2, 3, 4, 5]
    ring2 = [3, 4, 6, 7, 8]
    add_aromatic_ring(graph, atoms, ring1)
    add_aromatic_ring(graph, atoms, ring2)
    return [atom.id for atom in atoms]


def build_pyrene(graph: MolGraph) -> list[int]:
    template_path = os.path.join(
        os.path.dirname(__file__),
        "..",
        "src",
        "chemname",
        "templates",
        "fused",
        "pyrene_iupac2004.mol",
    )
    template = load_template(template_path)
    atoms = [graph.add_atom(atom.element, float(i), 0.0) for i, atom in enumerate(template.atoms)]
    for bond in template.bonds:
        graph.add_bond(
            atoms[bond.a].id,
            atoms[bond.b].id,
            order=1 if bond.aromatic else bond.order,
            is_aromatic=bond.aromatic,
        )
    return [atom.id for atom in atoms]


class ChemNamePR26Test(unittest.TestCase):
    def test_pyrene(self):
        graph = MolGraph()
        build_pyrene(graph)
        self.assertEqual(iupac_name(graph), "pyrene")

    def test_benzofuran(self):
        graph = MolGraph()
        build_benzofuran(graph)
        self.assertEqual(iupac_name(graph), "benzofuran")

    def test_chlorobenzofuran(self):
        graph = MolGraph()
        atoms = build_benzofuran(graph)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[7], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2-chlorobenzofuran")

    def test_benzothiophene(self):
        graph = MolGraph()
        build_benzothiophene(graph)
        self.assertEqual(iupac_name(graph), "benzothiophene")

    def test_chloropyrene(self):
        graph = MolGraph()
        atoms = build_pyrene(graph)
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_iupac2004.mol",
        )
        template = load_template(template_path)
        loc1_atom = next(idx for idx, loc in template.locant_by_atom_idx.items() if loc == 1)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[loc1_atom], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1-chloropyrene")

    def test_pyrene_dichloro_1_6(self):
        graph = MolGraph()
        atoms = build_pyrene(graph)
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_iupac2004.mol",
        )
        template = load_template(template_path)
        loc_to_idx = {loc: idx for idx, loc in template.locant_by_atom_idx.items()}
        for loc in (1, 6):
            cl = graph.add_atom("Cl", -1.0, 0.0)
            graph.add_bond(atoms[loc_to_idx[loc]], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1,6-dichloropyrene")

    def test_pyrene_dichloro_1_8(self):
        graph = MolGraph()
        atoms = build_pyrene(graph)
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_iupac2004.mol",
        )
        template = load_template(template_path)
        loc_to_idx = {loc: idx for idx, loc in template.locant_by_atom_idx.items()}
        for loc in (1, 8):
            cl = graph.add_atom("Cl", -1.0, 0.0)
            graph.add_bond(atoms[loc_to_idx[loc]], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "1,8-dichloropyrene")

    def test_pyrene_dichloro_2_7(self):
        graph = MolGraph()
        atoms = build_pyrene(graph)
        template_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_iupac2004.mol",
        )
        template = load_template(template_path)
        loc_to_idx = {loc: idx for idx, loc in template.locant_by_atom_idx.items()}
        for loc in (2, 7):
            cl = graph.add_atom("Cl", -1.0, 0.0)
            graph.add_bond(atoms[loc_to_idx[loc]], cl.id, order=1)
        self.assertEqual(iupac_name(graph), "2,7-dichloropyrene")

    def test_pyrene_scheme_switch(self):
        graph = MolGraph()
        atoms = build_pyrene(graph)
        template_path_iupac = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_iupac2004.mol",
        )
        template_path_cas = os.path.join(
            os.path.dirname(__file__),
            "..",
            "src",
            "chemname",
            "templates",
            "fused",
            "pyrene_cas.mol",
        )
        template_iupac = load_template(template_path_iupac)
        template_cas = load_template(template_path_cas)
        loc2_atom = next(idx for idx, loc in template_iupac.locant_by_atom_idx.items() if loc == 2)
        loc2_cas = template_cas.locant_by_atom_idx.get(loc2_atom)
        cl = graph.add_atom("Cl", -1.0, 0.0)
        graph.add_bond(atoms[loc2_atom], cl.id, order=1)
        name_iupac = iupac_name(graph, NameOptions(fused_numbering_scheme="iupac2004"))
        name_cas = iupac_name(graph, NameOptions(fused_numbering_scheme="cas"))
        self.assertEqual(name_iupac, "2-chloropyrene")
        self.assertNotEqual(name_iupac, name_cas)


if __name__ == "__main__":
    unittest.main()
