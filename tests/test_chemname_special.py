"""Pruebas de plantillas especiales (fase A: carbohidratos y esteroides)."""

import os
import sys
import unittest
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemio.rdkit_io import molfile_to_molgraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.chemname.molview import MolView
from chemuson.chemname.special import detect_special_template
from chemuson.utils.resources import open_resource_path


def _load_special_graph(filename: str) -> MolGraph:
    with open_resource_path("chemname", "templates", "special", filename) as path:
        molblock = Path(path).read_text(encoding="utf-8")
    return molfile_to_molgraph(molblock)


class ChemNameSpecialTemplatesTest(unittest.TestCase):
    """Cobertura MVP para plantillas de carbohidratos y esteroides."""

    def test_exact_carbohydrate_templates(self):
        self.assertEqual(
            iupac_name(_load_special_graph("alpha_d_glucopyranose.mol")),
            "alpha-d-glucopyranose",
        )
        self.assertEqual(
            iupac_name(_load_special_graph("beta_d_glucopyranose.mol")),
            "beta-d-glucopyranose",
        )
        self.assertEqual(
            iupac_name(_load_special_graph("beta_d_fructofuranose.mol")),
            "beta-d-fructofuranose",
        )
        self.assertEqual(
            iupac_name(_load_special_graph("d_ribose.mol")),
            "d-ribose",
        )

    def test_exact_steroid_templates(self):
        self.assertEqual(iupac_name(_load_special_graph("androstane_core.mol")), "androstane")
        self.assertEqual(iupac_name(_load_special_graph("cholestane_core.mol")), "cholestane")

    def test_special_substituted_hydroxy(self):
        graph = _load_special_graph("androstane_core.mol")
        anchor_id = min(graph.atoms.keys())
        hydroxyl = graph.add_atom("O", -2.0, 1.0, explicit_h=1)
        graph.add_bond(anchor_id, hydroxyl.id, order=1)
        name = iupac_name(graph)
        self.assertIn("hydroxyandrostane", name)

    def test_special_substituted_unknown_fallback(self):
        graph = _load_special_graph("androstane_core.mol")
        anchor_id = min(graph.atoms.keys())
        chlorine = graph.add_atom("Cl", -2.0, 1.0)
        graph.add_bond(anchor_id, chlorine.id, order=1)
        self.assertEqual(iupac_name(graph), "androstane (substituted)")

    def test_disable_special_templates_degrades(self):
        graph = _load_special_graph("alpha_d_glucopyranose.mol")
        name = iupac_name(
            graph,
            NameOptions(
                enable_special_templates=False,
                enable_experimental=True,
            ),
        )
        self.assertNotEqual(name, "alpha-d-glucopyranose")

    def test_exact_matching_works_without_rdkit(self):
        import chemuson.chemname.special as special_mod

        graph = _load_special_graph("d_ribose.mol")
        view = MolView(graph)
        previous = special_mod.Chem
        try:
            special_mod.Chem = None
            detected = detect_special_template(view)
        finally:
            special_mod.Chem = previous
        self.assertIsNotNone(detected)
        assert detected is not None
        self.assertEqual(detected[0], "d-ribose")

    def test_non_special_structure_keeps_standard_path(self):
        graph = MolGraph()
        c1 = graph.add_atom("C", 0.0, 0.0)
        c2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(c1.id, c2.id, order=1)
        self.assertEqual(iupac_name(graph), "ethane")


if __name__ == "__main__":
    unittest.main()

