"""Reference-backed aldehyde multiplicity and seniority checks."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    "smiles",
    ["O=CCC=O", "C(=O)CC=O"],
    ids=["registry-smiles", "reversed-atom-order"],
)
def test_terminal_dialdehydes_use_multiplicative_dial_suffix(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "propanedial"


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("O=CC(=O)O", "2-oxoethanoic acid"),
        ("O=CCC(=O)O", "3-oxopropanoic acid"),
        ("O=CCO", "2-hydroxyethanal"),
    ],
    ids=["acid-outranks-aldehyde", "acid-retains-priority", "aldehyde-outranks-alcohol"],
)
def test_dialdehyde_suffix_change_preserves_senior_group_priority(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected
