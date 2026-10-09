"""Reference-backed multiplicative alcohol suffix and seniority checks."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("OCCO", "ethane-1,2-diol"),
        ("OCCCO", "propane-1,3-diol"),
        ("OCCC(O)", "propane-1,3-diol"),
    ],
    ids=["ethane-diol", "propane-diol", "reversed-propane-diol"],
)
def test_multiple_alcohols_use_multiplicative_diol_suffix(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCO", "ethan-1-ol"),
        ("CC(=O)CCO", "4-hydroxybutan-2-one"),
        ("CC(O)C(=O)O", "2-hydroxypropanoic acid"),
    ],
    ids=["single-alcohol", "ketone-outranks-alcohol", "acid-outranks-alcohol"],
)
def test_diol_support_preserves_single_and_higher_priority_functions(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected
