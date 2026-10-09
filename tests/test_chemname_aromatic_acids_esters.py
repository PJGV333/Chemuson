"""Reference-backed benzoic-acid and benzoate naming cases."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("O=C(O)c1ccccc1", "benzoic acid"),
        ("COC(=O)c1ccccc1", "methyl benzoate"),
        ("CCOC(=O)c1ccccc1", "ethyl benzoate"),
    ],
    ids=["benzoic-acid", "methyl-benzoate", "ethyl-benzoate"],
)
def test_aromatic_carboxyl_names_use_reference_backed_forms(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


def test_phenylacetic_acid_is_not_misclassified_as_benzoic_acid() -> None:
    graph = smiles_to_molgraph("O=C(O)Cc1ccccc1")
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "2-phenylethanoic acid"


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        ("O=C(O)c1ccccc1", "c1ccccc1C(=O)O", "benzoic acid"),
        ("COC(=O)c1ccccc1", "c1ccccc1C(=O)OC", "methyl benzoate"),
    ],
    ids=["benzoic-acid-order", "benzoate-order"],
)
def test_aromatic_carboxyl_names_are_invariant_to_atom_order(
    first: str, second: str, expected: str
) -> None:
    options = NameOptions(rdkit_isolated=False)
    assert iupac_name(smiles_to_molgraph(first), options) == expected
    assert iupac_name(smiles_to_molgraph(second), options) == expected
