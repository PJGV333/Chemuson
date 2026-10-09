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


def test_reference_corpus_has_78_independently_resolved_structures() -> None:
    assert len(CASES) == 78
    assert len({case["pubchem_cid"] for case in CASES}) == 78
    assert all(case["pubchem_name"] and case["pubchem_connectivity_smiles"] for case in CASES)
    assert all(case["formula"] and case["smiles"] for case in CASES)
    assert {case["baseline_classification"] for case in CASES} == {
        "correct",
        "incorrect",
        "unsupported",
        "reference_pending",
    }
    assert all(
        bool(case["target_name"]) == (case["baseline_classification"] != "reference_pending")
        for case in CASES
    )


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
