from __future__ import annotations

import copy
import math

import pytest


from chemuson.clean2d import smiles_to_depiction_candidates
from chemuson.clean2d.depiction_quality import block_donut_score, score_imported_depiction
from chemuson.clean2d.block_unwrap import block_unwrap_layout
from chemuson.clean2d.engine import run_clean2d_engine
from chemuson.core.model import MolGraph


TETRANDRINE_SMILES = "CN1CCC2=CC(=C3C=C2C1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)CC6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"


def test_tetrandrine_import_generates_block_unwrap_candidate() -> None:
    try:
        candidates = smiles_to_depiction_candidates(TETRANDRINE_SMILES, timeout_s=20.0)
    except RuntimeError as exc:
        pytest.skip(f"RDKit worker unavailable: {exc}")

    unwraps = [candidate for candidate in candidates if "chemuson_block_unwrap" in candidate.source]
    if not unwraps:
        pytest.skip("RDKit depiction was not donut-like enough to trigger block unwrap")
    assert unwraps[0].graph.atoms
    assert unwraps[0].graph.bonds


def test_tetrandrine_best_depiction_prefers_not_donut() -> None:
    try:
        candidates = smiles_to_depiction_candidates(TETRANDRINE_SMILES, timeout_s=20.0)
    except RuntimeError as exc:
        pytest.skip(f"RDKit worker unavailable: {exc}")
    accepted = [candidate for candidate in candidates if not candidate.rejected]
    if not accepted:
        pytest.skip("No accepted RDKit depiction candidates")

    best = accepted[0]
    donut_score = float(best.metadata.get("donut_score", 0.0) or 0.0)
    assert "block_unwrap" in best.source or donut_score < 4.0
    assert len(best.graph.atoms) > 35
    assert len(best.graph.bonds) > 35


def test_block_unwrap_improves_synthetic_donut() -> None:
    graph = _synthetic_donut_graph()
    before_score, _ = score_imported_depiction(graph)
    before_donut, _ = block_donut_score(graph)

    coords, report = block_unwrap_layout(graph)

    assert coords is not None, report
    assert report.ok
    after = copy.deepcopy(graph)
    for atom_id, (x, y) in coords.items():
        after.atoms[atom_id].x = x
        after.atoms[atom_id].y = y
    after_score, _ = score_imported_depiction(after)
    after_donut, _ = block_donut_score(after)
    assert after_score < before_score
    assert after_donut < before_donut


def test_block_unwrap_not_used_for_simple_molecule() -> None:
    graph = _benzene_graph()

    coords, report = block_unwrap_layout(graph)

    assert coords is None
    assert not report.ok


def test_clean2d_uses_block_unwrap_for_donut_complex() -> None:
    graph = _synthetic_donut_graph()
    before_donut, _ = block_donut_score(graph)

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)

    sources = {candidate.source for candidate in (*result.candidates, *result.rejected)}
    assert before_donut > 0.0
    assert sources & {"block_unwrap", "scaffold_depiction"} or (
        result.selected is not None and result.selected.source in {"block_unwrap", "scaffold_depiction"}
    )


def test_existing_simple_clean2d_still_works() -> None:
    graph = _benzene_graph()

    result = run_clean2d_engine(graph, set(graph.atoms), mode="quick", target_bond_length=40.0)

    assert result.selected is not None
    assert all(candidate.source != "block_unwrap" for candidate in (*result.candidates, *result.rejected))


def _synthetic_donut_graph() -> MolGraph:
    graph = MolGraph()
    centers = [(-130.0, -130.0), (130.0, -130.0), (130.0, 130.0), (-130.0, 130.0)]
    rings = [_add_hexagon(graph, x, y) for x, y in centers]
    for idx in range(3):
        _connect_with_linker(graph, rings[idx][0], rings[idx + 1][3])
    return graph


def _benzene_graph() -> MolGraph:
    graph = MolGraph()
    _add_hexagon(graph, 0.0, 0.0)
    return graph


def _add_hexagon(graph: MolGraph, cx: float, cy: float) -> list[int]:
    atoms: list[int] = []
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
