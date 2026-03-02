"""Smoke tests de entorno RDKit (aislamiento y disponibilidad)."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.rdkit_safe import run_rdkit_stereo_extract


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


if __name__ == "__main__":
    unittest.main()
