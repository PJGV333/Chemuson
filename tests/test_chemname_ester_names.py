"""Reference-backed ester naming: both ester components must be represented."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.chemname.errors import ChemNameNotSupported


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("COC(=O)C", "methyl acetate"),
        ("CCOC(=O)C", "ethyl acetate"),
        ("CCC(=O)OC", "methyl propanoate"),
        ("CCOC(=O)CC", "ethyl propanoate"),
        ("CC(C)C(=O)OC", "methyl 2-methylpropanoate"),
        ("C=CC(=O)OC", "methyl prop-2-enoate"),
        ("CCOC(=O)C(O)C", "ethyl 2-hydroxypropanoate"),
        ("CCOC(=O)CCN", "ethyl 3-aminopropanoate"),
        ("COC(=O)C(N)C(=O)C", "methyl 2-amino-3-oxobutanoate"),
    ],
    ids=[
        "methyl-acetate",
        "ethyl-acetate",
        "methyl-propanoate",
        "ethyl-propanoate",
        "branched-acid-ester",
        "unsaturated-acid-ester",
        "hydroxy-ester",
        "amino-ester",
        "amino-oxo-ester",
    ],
)
def test_ester_name_retains_organyl_and_acid_components(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(=O)O", "ethanoic acid"),
        ("CC(=O)CC(=O)O", "3-oxobutanoic acid"),
    ],
    ids=["carboxylic-acid-is-not-an-ester", "acid-still-outranks-ketone"],
)
def test_ester_component_does_not_override_nonester_functional_groups(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "[13CH3]OC(=O)C",
        "CCC[C@H](C)OC(=O)C",
        "CCOCCOC(=O)C",
    ],
    ids=["isotopic-organyl", "stereogenic-organyl", "heteroatom-decorated-organyl"],
)
def test_unsupported_ester_organyl_information_fails_closed(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "N/D"
    with pytest.raises(ChemNameNotSupported):
        iupac_name(
            graph,
            NameOptions(rdkit_isolated=False, return_nd_on_fail=False),
        )


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        ("CCOC(=O)C", "CC(=O)OCC", "ethyl acetate"),
        (
            "COC(=O)C(N)C(=O)C",
            "CC(=O)C(N)C(=O)OC",
            "methyl 2-amino-3-oxobutanoate",
        ),
    ],
    ids=["simple-ester", "multifunctional-ester"],
)
def test_ester_naming_is_invariant_to_smiles_atom_order(
    first: str, second: str, expected: str
) -> None:
    options = NameOptions(rdkit_isolated=False)
    assert iupac_name(smiles_to_molgraph(first), options) == expected
    assert iupac_name(smiles_to_molgraph(second), options) == expected
