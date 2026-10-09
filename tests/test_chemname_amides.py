"""Reference-backed N-substituted and multiple amide regressions."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.chemname.errors import ChemNameNotSupported


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CNC(C)=O", "N-methylethanamide"),
        ("CCNC(C)=O", "N-ethylethanamide"),
    ],
    ids=["n-methyl-systematic", "n-ethyl-systematic"],
)
def test_amide_name_retains_n_alkyl_substituent(smiles: str, expected: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


def test_two_principal_amide_groups_use_diamide_suffix() -> None:
    graph = smiles_to_molgraph("NC(=O)CC(=O)N")
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "propanediamide"


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(=O)N", "ethanamide"),
        ("CC(=O)CC(=O)N", "3-oxobutanamide"),
    ],
    ids=["primary-amide-control", "amide-outranks-ketone"],
)
def test_n_substitution_support_preserves_primary_and_mixed_amides(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    "smiles",
    ["[13CH3]NC(C)=O", "CCC[C@H](C)NC(C)=O"],
    ids=["isotopic-n-alkyl", "stereogenic-n-alkyl"],
)
def test_unsupported_n_alkyl_metadata_fails_closed(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "N/D"
    with pytest.raises(ChemNameNotSupported):
        iupac_name(
            graph,
            NameOptions(rdkit_isolated=False, return_nd_on_fail=False),
        )
