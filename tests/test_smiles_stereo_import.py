from __future__ import annotations

import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.rdkit_io import smiles_to_molgraph_best_depiction
from chemuson.core.model import BondStyle


TETRANDRINE_SMILES = "CN1CCC2=CC(=C3C=C2C1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)CC6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"
VANCOMYCIN_SMILES = "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"


def test_chiral_smiles_import_creates_wedge_or_hash() -> None:
    graph = _import_or_skip("C[C@H](O)F")
    assert _wedge_hash_count(graph) >= 1
    assert any(atom.stereo_cip for atom in graph.atoms.values())


def test_amino_acid_chiral_smiles_import_creates_wedge_or_hash() -> None:
    graph = _import_or_skip("N[C@@H](C)C(=O)O")
    assert _wedge_hash_count(graph) >= 1
    assert any(atom.stereo_cip for atom in graph.atoms.values())


def test_tetrandrine_import_preserves_visual_stereo() -> None:
    graph = _import_or_skip(TETRANDRINE_SMILES, timeout_s=20.0)
    assert _wedge_hash_count(graph) >= 1


def test_vancomycin_import_preserves_visual_stereo() -> None:
    graph = _import_or_skip(VANCOMYCIN_SMILES, timeout_s=20.0)
    assert _wedge_hash_count(graph) >= 4


def _import_or_skip(smiles: str, timeout_s: float = 10.0):
    try:
        return smiles_to_molgraph_best_depiction(smiles, timeout_s=timeout_s)
    except RuntimeError as exc:
        pytest.skip(f"RDKit worker unavailable: {exc}")


def _wedge_hash_count(graph) -> int:
    return sum(1 for bond in graph.bonds.values() if bond.style in {BondStyle.WEDGE, BondStyle.HASHED})
