"""Reference-backed ketone multiplicity and suffix-seniority checks."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    "smiles",
    ["CC(=O)CCC(=O)C", "O=C(C)CCC(=O)C"],
    ids=["registry-smiles", "alternate-atom-order"],
)
def test_two_principal_ketones_use_dione_suffix(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "hexane-2,5-dione"


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(=O)CC", "butan-2-one"),
        ("CC(=O)CC(=O)O", "3-oxobutanoic acid"),
    ],
    ids=["single-ketone", "carboxylic-acid-outranks-ketone"],
)
def test_dione_support_does_not_degrade_single_or_mixed_functions(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected
