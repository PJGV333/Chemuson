"""Pruebas unitarias para test_chemname_pr0."""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemname import iupac_name


class ChemNamePR0Test(unittest.TestCase):
    """Casos de prueba para ChemNamePR0Test."""
    def test_stub_returns_nd(self):
        """Verifica stub returns nd.

        Returns:
            None.

        """
        self.assertEqual(iupac_name(None), "N/D")


if __name__ == "__main__":
    unittest.main()
