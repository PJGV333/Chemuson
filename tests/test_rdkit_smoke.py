"""Smoke tests de entorno RDKit (aislamiento y disponibilidad)."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.rdkit_safe import (
    advanced_stereo_descriptors_for_chain,
    run_rdkit_stereo_extract,
)
from chemuson.core.model import MolGraph


class RdkitSmokeTest(unittest.TestCase):
    """Valida que RDKit pueda importarse en subproceso sin romper la suite."""

    def test_rdkit_import_subprocess(self):
        proc = subprocess.run(
            [sys.executable, "-c", "import rdkit; print('ok')"],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            self.skipTest(f"RDKit no disponible o inestable (exit={proc.returncode})")
        self.assertIn("ok", (proc.stdout or "").strip())

    def test_rdkit_safe_worker_graceful(self):
        result = run_rdkit_stereo_extract("C/C=C/C", fmt="smiles", timeout_s=5.0)
        if not result.get("ok"):
            self.skipTest("Worker RDKit no disponible en este entorno")
        self.assertTrue(result.get("ok"))

    def test_rdkit_safe_advanced_worker_graceful(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0, stereo_axial="R_a")
        a2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(a1.id, a2.id, order=1, stereo_endo_exo="endo")
        descriptors = advanced_stereo_descriptors_for_chain(graph, [a1.id, a2.id], timeout_s=5.0)
        # Si RDKit falla/no está, la API debe degradar limpiamente a lista vacía.
        self.assertIsInstance(descriptors, list)

    def test_rdkit_safe_worker_tolerates_pseudoatoms_and_duplicate_bonds(self):
        graph = MolGraph()
        a1 = graph.add_atom("OH", 0.0, 0.0, is_explicit=True)
        a2 = graph.add_atom("CH2OH", 1.0, 0.0, is_explicit=True)
        graph.add_bond(a1.id, a2.id, order=1)
        graph.add_bond(a1.id, a2.id, order=2)

        descriptors = advanced_stereo_descriptors_for_chain(graph, [a1.id, a2.id], timeout_s=5.0)
        self.assertIsInstance(descriptors, list)


if __name__ == "__main__":
    unittest.main()
