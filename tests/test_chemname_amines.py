"""Reference-backed diamine suffix and amine-seniority checks."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    "smiles",
    ["NCCN", "C(CN)N"],
    ids=["registry-smiles", "alternate-atom-order"],
)
def test_two_principal_amines_use_multiplicative_diamine_suffix(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "ethane-1,2-diamine"


def test_diamine_support_preserves_acid_suffix_seniority() -> None:
    graph = smiles_to_molgraph("CC(N)C(=O)O")
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "2-aminopropanoic acid"
