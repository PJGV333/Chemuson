"""Reference-backed simple aryl ketones and explicit unsupported boundaries."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC(=O)c1ccccc1", "1-phenylethan-1-one"),
        ("CC(=O)c1ccc(O)cc1", "1-(4-hydroxyphenyl)ethan-1-one"),
    ],
    ids=["acetophenone", "p-hydroxyacetophenone"],
)
def test_simple_aryl_ketones_use_reference_backed_names(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        ("CC(=O)c1ccccc1", "c1ccccc1C(=O)C", "1-phenylethan-1-one"),
        (
            "CC(=O)c1ccc(O)cc1",
            "Oc1ccc(C(C)=O)cc1",
            "1-(4-hydroxyphenyl)ethan-1-one",
        ),
    ],
    ids=["acetophenone-order", "hydroxyketone-order"],
)
def test_aryl_ketone_names_are_invariant_to_atom_order(
    first: str, second: str, expected: str
) -> None:
    options = NameOptions(rdkit_isolated=False)
    assert iupac_name(smiles_to_molgraph(first), options) == expected
    assert iupac_name(smiles_to_molgraph(second), options) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)c1cccc2ccccc12",
        "CC(=O)c1ccc([N+](C)(C)C)cc1",
        "CC(=O)c1ccc([13CH3])cc1",
        "CC(=O)c1ccc(C[C@H](O)C)cc1",
        "CC(=O)c1cc(N)ccc1C",
    ],
    ids=[
        "fused-ring",
        "charged-decoration",
        "isotopic-decoration",
        "stereogenic-decoration",
        "other-direct-aryl-substituents-out-of-scope",
    ],
)
def test_unsupported_aryl_ketone_decorations_fail_closed(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == "N/D"
