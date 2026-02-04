import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemname import iupac_name


class ChemNamePR0Test(unittest.TestCase):
    def test_stub_returns_nd(self):
        self.assertEqual(iupac_name(None), "N/D")


if __name__ == "__main__":
    unittest.main()
