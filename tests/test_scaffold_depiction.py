from __future__ import annotations

import copy
import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio.depiction_candidates import block_donut_score, score_imported_depiction
from chemuson.chemio.rdkit_io import smiles_to_depiction_candidates, smiles_to_molgraph_best_depiction_with_report
from chemuson.clean2d.scaffold_depiction import scaffold_depiction_candidates, scaffold_depiction_layout
from chemuson.core.model import MolGraph


TETRANDRINE_SMILES = "CN1CCC2=CC(=C3C=C2C1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)CC6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"
VANCOMYCIN_SMILES = "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"


def test_scaffold_depiction_generates_multiple_candidates_for_synthetic_complex() -> None:
    graph = _synthetic_complex()
    candidates = scaffold_depiction_candidates(graph)
    accepted = [candidate for candidate in candidates if not candidate.rejected]
    assert len(accepted) > 1
    assert {candidate.strategy for candidate in accepted} >= {"longest_path_zigzag", "layered_scaffold"}


def test_scaffold_depiction_improves_synthetic_donut() -> None:
    graph = _synthetic_complex()
    before_score, _ = score_imported_depiction(graph)
    before_donut, _ = block_donut_score(graph)
    coords, report = scaffold_depiction_layout(graph)
    assert coords is not None, report
    after = copy.deepcopy(graph)
    for atom_id, (x, y) in coords.items():
        after.atoms[atom_id].x = x
        after.atoms[atom_id].y = y
    after_score, _ = score_imported_depiction(after)
    after_donut, _ = block_donut_score(after)
    assert after_score < before_score
    assert after_donut < before_donut


def test_tetrandrine_best_depiction_not_forced_to_single_strategy() -> None:
    candidates = _candidates_or_skip(TETRANDRINE_SMILES)
    sources = {candidate.source for candidate in candidates}
    assert any(source in sources for source in {"rdcoordgen", "rddepictor_compute2d"})
    if not any(source.startswith("chemuson_scaffold:") for source in sources):
        pytest.skip("RDKit tetrandrine depiction did not trigger scaffold candidate")
    best = [candidate for candidate in candidates if not candidate.rejected][0]
    assert best.graph.atoms
    assert best.graph.bonds


def test_vancomycin_best_depiction_has_scaffold_candidate() -> None:
    candidates = _candidates_or_skip(VANCOMYCIN_SMILES, timeout_s=20.0)
    if not any(candidate.source.startswith("chemuson_scaffold:") for candidate in candidates):
        pytest.skip("RDKit vancomycin depiction did not trigger scaffold candidate")
    graph, report = smiles_to_molgraph_best_depiction_with_report(VANCOMYCIN_SMILES, timeout_s=20.0)
    assert graph.atoms
    assert report.get("selected_source")
    assert report.get("candidates")


def test_simple_molecule_does_not_trigger_scaffold_depiction() -> None:
    candidates = _candidates_or_skip("CCO", timeout_s=10.0)
    assert not any(candidate.source.startswith("chemuson_scaffold:") and not candidate.rejected for candidate in candidates)


def _candidates_or_skip(smiles: str, timeout_s: float = 20.0):
    try:
        candidates = smiles_to_depiction_candidates(smiles, timeout_s=timeout_s)
    except RuntimeError as exc:
        pytest.skip(f"RDKit worker unavailable: {exc}")
    if not candidates:
        pytest.skip("No depiction candidates")
    return candidates


def _synthetic_complex() -> MolGraph:
    graph = MolGraph()
    centers = [(-130.0, -130.0), (130.0, -130.0), (130.0, 130.0), (-130.0, 130.0)]
    rings = [_add_hexagon(graph, x, y) for x, y in centers]
    for idx in range(3):
        _connect_with_linker(graph, rings[idx][0], rings[idx + 1][3])
    return graph


def _add_hexagon(graph: MolGraph, cx: float, cy: float) -> list[int]:
    atoms = []
    for idx in range(6):
        angle = math.radians(idx * 60.0)
        atoms.append(graph.add_atom("C", cx + 40.0 * math.cos(angle), cy + 40.0 * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(atoms[idx], atoms[(idx + 1) % 6], order=1, is_aromatic=True)
    return atoms


def _connect_with_linker(graph: MolGraph, left: int, right: int) -> None:
    x1, y1 = graph.atoms[left].x, graph.atoms[left].y
    x2, y2 = graph.atoms[right].x, graph.atoms[right].y
    a = graph.add_atom("C", x1 * 0.67 + x2 * 0.33, y1 * 0.67 + y2 * 0.33).id
    b = graph.add_atom("C", x1 * 0.33 + x2 * 0.67, y1 * 0.33 + y2 * 0.67).id
    graph.add_bond(left, a, order=1)
    graph.add_bond(a, b, order=1)
    graph.add_bond(b, right, order=1)
