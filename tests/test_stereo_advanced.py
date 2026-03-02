"""Pruebas de estereoquímica avanzada (fase B, best-effort)."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.chemname.molview import MolView
from chemuson.chemname.render import render_name
from chemuson.chemname.stereo_advanced import extract_advanced_stereo

try:
    from chemuson.chemio.rdkit_io import smiles_to_molgraph
    from rdkit import Chem  # noqa: F401

    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False


def _build_linear_chain(graph: MolGraph, length: int) -> list[int]:
    atoms = [graph.add_atom("C", float(i), 0.0) for i in range(length)]
    for i in range(length - 1):
        graph.add_bond(atoms[i].id, atoms[i + 1].id, order=1)
    return [atom.id for atom in atoms]


class StereoAdvancedTest(unittest.TestCase):
    """Cobertura mínima de axialidad/helicidad/endo-exo/si-re."""

    def test_render_priority_helical_then_axial_then_classic(self):
        name = render_name(
            [],
            "butane",
            stereo_descriptors=["2R", "2E", "1R_a", "M"],
        )
        self.assertEqual(name, "(M,1R_a,2E,2R)-butane")

    def test_extract_annotations_for_advanced_descriptors(self):
        graph = MolGraph()
        chain = _build_linear_chain(graph, 4)
        graph.get_atom(chain[1]).stereo_axial = "Ra"
        graph.get_atom(chain[2]).stereo_si_re = "Si"
        graph.get_atom(chain[3]).stereo_helical = "P"
        graph.find_bond_between(chain[1], chain[2]).stereo_endo_exo = "exo"

        descriptors = extract_advanced_stereo(
            MolView(graph),
            chain=chain,
            opts=NameOptions(rdkit_isolated=False),
        )
        self.assertIn("P", descriptors)
        self.assertIn("2R_a", descriptors)
        self.assertIn("3si", descriptors)
        self.assertIn("2exo", descriptors)

    def test_iupac_uses_advanced_descriptors_when_enabled(self):
        graph = MolGraph()
        chain = _build_linear_chain(graph, 3)
        graph.get_atom(chain[1]).stereo_axial = "S_a"
        name = iupac_name(
            graph,
            NameOptions(
                enable_experimental=True,
                enable_advanced_stereo=True,
                rdkit_isolated=False,
            ),
        )
        self.assertIn("2S_a", name)

    def test_iupac_omits_advanced_descriptors_when_disabled(self):
        graph = MolGraph()
        chain = _build_linear_chain(graph, 3)
        graph.get_atom(chain[1]).stereo_axial = "S_a"
        name = iupac_name(
            graph,
            NameOptions(
                enable_experimental=True,
                enable_advanced_stereo=False,
                rdkit_isolated=False,
            ),
        )
        self.assertNotIn("S_a", name)

    def test_rdkit_isolation_failure_degrades_without_exception(self):
        import chemuson.chemname.stereo_advanced as stereo_mod

        graph = MolGraph()
        chain = _build_linear_chain(graph, 3)
        graph.get_atom(chain[1]).stereo_axial = "R_a"
        previous = stereo_mod.safe_advanced_stereo_descriptors_for_chain
        try:
            def _boom(*_args, **_kwargs):
                raise RuntimeError("worker failed")

            stereo_mod.safe_advanced_stereo_descriptors_for_chain = _boom
            descriptors = extract_advanced_stereo(
                MolView(graph),
                chain=chain,
                opts=NameOptions(rdkit_isolated=True),
            )
        finally:
            stereo_mod.safe_advanced_stereo_descriptors_for_chain = previous
        self.assertIn("2R_a", descriptors)

    @unittest.skipIf(not RDKIT_AVAILABLE, "RDKit no disponible")
    def test_allene_axial_best_effort(self):
        graph = smiles_to_molgraph("F/C=C=C/F")
        name = iupac_name(
            graph,
            NameOptions(
                enable_experimental=True,
                enable_advanced_stereo=True,
                rdkit_isolated=True,
            ),
        )
        self.assertIsInstance(name, str)
        self.assertNotEqual(name, "")

    @unittest.skipIf(not RDKIT_AVAILABLE, "RDKit no disponible")
    def test_biaryl_atropisomeric_best_effort(self):
        graph = smiles_to_molgraph("c1ccc(cc1)-c1ccccc1")
        name = iupac_name(
            graph,
            NameOptions(
                enable_experimental=True,
                enable_advanced_stereo=True,
                rdkit_isolated=True,
            ),
        )
        self.assertIsInstance(name, str)
        self.assertNotEqual(name, "")


if __name__ == "__main__":
    unittest.main()

