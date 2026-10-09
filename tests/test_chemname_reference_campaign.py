"""Independent PubChem-backed reference-corpus integrity and baseline tests."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "tests/data/chemname_iupac_reference_campaign.psv"


def _load_cases() -> list[dict[str, str]]:
    with CORPUS.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="|"))


CASES = _load_cases()


def test_reference_corpus_has_80_independently_resolved_structures() -> None:
    assert len(CASES) == 80
    assert len({case["pubchem_cid"] for case in CASES}) == 80
    assert all(case["pubchem_name"] and case["pubchem_connectivity_smiles"] for case in CASES)
    assert all(case["formula"] and case["smiles"] for case in CASES)
    assert {case["baseline_classification"] for case in CASES} == {
        "correct",
        "incorrect",
        "unsupported",
        "reference_pending",
    }
    assert all(case["target_name"] for case in CASES)
    assert [case["id"] for case in CASES if case["baseline_classification"] == "reference_pending"] == [
        "ester_ethyl_ethanoate"
    ]


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_pubchem_connectivity_reference_matches_exact_input(case: dict[str, str]) -> None:
    rdkit = pytest.importorskip("rdkit.Chem")
    input_mol = rdkit.MolFromSmiles(case["smiles"])
    reference_mol = rdkit.MolFromSmiles(case["pubchem_connectivity_smiles"])
    assert input_mol is not None
    assert reference_mol is not None
    assert rdkit.MolToSmiles(input_mol, canonical=True) == rdkit.MolToSmiles(
        reference_mol, canonical=True
    )


@pytest.mark.parametrize(
    "case",
    [case for case in CASES if case["baseline_classification"] == "correct"],
    ids=[case["id"] for case in CASES if case["baseline_classification"] == "correct"],
)
def test_reference_backed_correct_names_remain_stable(case: dict[str, str]) -> None:
    graph = smiles_to_molgraph(case["smiles"])
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == case["baseline_name"]


SYSTEMATIC_VARIANTS = {
    "amide_nmethylethanamide": "N-methylethanamide",
    "amide_nethylethanamide": "N-ethylethanamide",
}


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_final_campaign_names_match_targets_or_adjudicated_systematic_variants(
    case: dict[str, str],
) -> None:
    graph = smiles_to_molgraph(case["smiles"])
    actual = iupac_name(graph, NameOptions(rdkit_isolated=False))
    expected = SYSTEMATIC_VARIANTS.get(case["id"], case["target_name"])
    assert actual == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("O=C(O)C(=O)O", "ethanedioic acid"),
        ("O=C(O)CCC(=O)O", "butanedioic acid"),
    ],
    ids=["ethanedioic-acid", "butanedioic-acid"],
)
def test_all_principal_carboxyl_groups_use_multiplicative_suffix(
    smiles: str, expected: str
) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph, NameOptions(rdkit_isolated=False)) == expected
