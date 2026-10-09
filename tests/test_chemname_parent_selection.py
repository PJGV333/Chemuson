"""Reference-backed parent-chain selection and seniority regressions."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(C)C(=O)O", "2-methylpropanoic acid"),
        ("CC(C)C=O", "2-methylpropanal"),
    ],
    ids=["branched-carboxylic-acid", "branched-aldehyde"],
)
def test_principal_carbonyl_group_is_included_in_selected_parent(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(=O)CC(=O)O", "3-oxobutanoic acid"),
        ("O=CCO", "2-hydroxyethanal"),
    ],
    ids=["acid-outranks-ketone", "aldehyde-outranks-alcohol"],
)
def test_parent_selection_keeps_functional_group_seniority(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        ("CC(C)C(=O)O", "O=C(O)C(C)C", "2-methylpropanoic acid"),
        ("CC(C)C=O", "O=CC(C)C", "2-methylpropanal"),
    ],
    ids=["acid-atom-order", "aldehyde-atom-order"],
)
def test_parent_selection_is_invariant_to_smiles_atom_order(
    first: str, second: str, expected: str
) -> None:
    options = NameOptions(rdkit_isolated=False)
    first_name = iupac_name(smiles_to_molgraph(first), options)
    second_name = iupac_name(smiles_to_molgraph(second), options)
    assert first_name == expected
    assert second_name == expected
