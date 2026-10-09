"""Atom-order stability for mixed and symmetric aromatic locant choices."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        (
            "Cc1ccc(O)cc1",
            "Oc1ccc(C)cc1",
            "1-hydroxy-4-methylbenzene",
        ),
        (
            "Oc1ccc([N+](=O)[O-])cc1",
            "[O-][N+](=O)c1ccc(O)cc1",
            "1-hydroxy-4-nitrobenzene",
        ),
        (
            "Cc1cc(C)cc(C)c1",
            "c1c(C)cc(C)cc1C",
            "1,3,5-trimethylbenzene",
        ),
    ],
    ids=["mixed-hydroxy-methyl", "mixed-hydroxy-nitro", "symmetric-trisubstitution"],
)
def test_aromatic_locant_ties_are_invariant_to_atom_order(
    first: str, second: str, expected: str
) -> None:
    options = NameOptions(rdkit_isolated=False)
    assert iupac_name(smiles_to_molgraph(first), options) == expected
    assert iupac_name(smiles_to_molgraph(second), options) == expected
