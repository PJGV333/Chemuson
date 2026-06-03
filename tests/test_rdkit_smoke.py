"""Smoke tests de entorno RDKit (aislamiento y disponibilidad)."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.rdkit_safe import (
    advanced_stereo_descriptors_for_chain,
    molgraph_to_smiles_isolated,
    smiles_to_molgraph_isolated,
    run_rdkit_stereo_extract,
    text_to_molblock,
)
import chemuson.chemio.rdkit_safe as rdkit_safe
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

    def test_rdkit_worker_probe_does_not_import_rdkit_in_parent(self):
        env = os.environ.copy()
        src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
        env["PYTHONPATH"] = src_path + os.pathsep + env.get("PYTHONPATH", "")
        script = (
            "import sys\n"
            "from chemuson.chemio.rdkit_safe import is_rdkit_worker_available\n"
            "before = any(name == 'rdkit' or name.startswith('rdkit.') for name in sys.modules)\n"
            "is_rdkit_worker_available(timeout_s=5.0)\n"
            "after = any(name == 'rdkit' or name.startswith('rdkit.') for name in sys.modules)\n"
            "print(f'{before} {after}')\n"
            "raise SystemExit(0 if (not before and not after) else 1)\n"
        )
        proc = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("False False", proc.stdout)

    def test_rdkit_safe_advanced_worker_graceful(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0, stereo_axial="R_a")
        a2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(a1.id, a2.id, order=1, stereo_endo_exo="endo")
        descriptors = advanced_stereo_descriptors_for_chain(graph, [a1.id, a2.id], timeout_s=5.0)
        # Si RDKit falla/no está, la API debe degradar limpiamente a lista vacía.
        self.assertIsInstance(descriptors, list)

    def test_rdkit_safe_text_to_molblock_and_graph(self):
        payload = text_to_molblock("CC", fmt="smiles", timeout_s=5.0)
        if not payload.get("ok"):
            self.skipTest("Worker RDKit no disponible en este entorno")
        self.assertIn("molblock", payload)
        graph, error = smiles_to_molgraph_isolated("CC", timeout_s=5.0)
        self.assertIsNone(error)
        self.assertIsNotNone(graph)

    def test_rdkit_safe_molgraph_to_smiles_payload_uses_chain_atom_ids(self):
        graph = MolGraph()
        a1 = graph.add_atom("C", 0.0, 0.0)
        a2 = graph.add_atom("C", 1.0, 0.0)
        graph.add_bond(a1.id, a2.id, order=1)
        captured = {}
        original_run_worker = rdkit_safe._run_worker

        def _fake_run_worker(request, timeout_s):
            captured["request"] = request
            captured["timeout_s"] = timeout_s
            return {"ok": True, "smiles": "CC"}

        rdkit_safe._run_worker = _fake_run_worker
        try:
            smiles, error = molgraph_to_smiles_isolated(graph, timeout_s=1.5)
        finally:
            rdkit_safe._run_worker = original_run_worker

        self.assertEqual(smiles, "CC")
        self.assertIsNone(error)
        self.assertEqual(captured["request"]["mode"], "graph_to_smiles")
        self.assertEqual(captured["request"]["chain"], [])
        self.assertEqual(captured["timeout_s"], 1.5)

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
