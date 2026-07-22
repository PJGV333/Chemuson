from __future__ import annotations

import math

import pytest


from chemuson.chemio import rdkit_io, rdkit_safe
from chemuson.chemio.depiction_candidates import DepictionCandidate, score_imported_depiction
from chemuson.core.model import MolGraph


VANCOMYCIN_SMILES = "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"


def _simple_graph(offset: float) -> MolGraph:
    graph = MolGraph()
    a1 = graph.add_atom("C", offset, 0.0)
    a2 = graph.add_atom("C", offset + 40.0, 0.0)
    graph.add_bond(a1.id, a2.id, order=1)
    return graph


def test_smiles_depiction_candidates_uses_worker_payload(monkeypatch) -> None:
    captured = {}

    def fake_worker(request, timeout_s):
        captured["request"] = request
        captured["timeout_s"] = timeout_s
        return {
            "ok": True,
            "candidates": [
                {"source": "slow", "ok": True, "molblock": "mol-a", "metadata": {}},
                {"source": "fast", "ok": True, "molblock": "mol-b", "metadata": {}},
            ],
        }

    def fake_parser(molblock: str) -> MolGraph:
        return _simple_graph(0.0 if molblock == "mol-a" else 100.0)

    def fake_score(graph: MolGraph, *, target_bond_length: float = 40.0):
        first_x = min(atom.x for atom in graph.atoms.values())
        return (20.0 if first_x == 0.0 else 5.0), {"first_x": first_x}

    monkeypatch.setattr(rdkit_safe, "_run_worker", fake_worker)
    monkeypatch.setattr(rdkit_io, "molfile_to_molgraph", fake_parser)
    monkeypatch.setattr(rdkit_io, "score_imported_depiction", fake_score)

    candidates = rdkit_io.smiles_to_depiction_candidates("CC", timeout_s=3.0)

    assert captured["request"]["mode"] == "smiles_depict_candidates"
    assert captured["request"]["smiles"] == "CC"
    assert captured["timeout_s"] == 3.0
    assert [candidate.source for candidate in candidates] == ["fast", "slow"]


def test_smiles_to_molgraph_best_depiction_prefers_lower_score(monkeypatch) -> None:
    best = _simple_graph(200.0)
    worse = _simple_graph(0.0)
    monkeypatch.setattr(
        rdkit_io,
        "smiles_to_depiction_candidates",
        lambda *_args, **_kwargs: (
            DepictionCandidate("worse", worse, 10.0),
            DepictionCandidate("best", best, 1.0),
        ),
    )

    graph = rdkit_io.smiles_to_molgraph_best_depiction("CC")

    assert graph is best


def test_depiction_candidate_contract_and_sorting_preserve_rejected_last(monkeypatch) -> None:
    accepted = _simple_graph(0.0)
    monkeypatch.setattr(
        rdkit_io,
        "smiles_to_depiction_candidates",
        lambda *_args, **_kwargs: (
            DepictionCandidate("zeta", accepted, 2.0, metadata={"origin": "worker"}),
            DepictionCandidate("alpha", MolGraph(), math.inf, True, "empty_molblock"),
            DepictionCandidate("beta", _simple_graph(100.0), 1.0),
        ),
    )

    graph, report = rdkit_io.smiles_to_molgraph_best_depiction_with_report("CC")

    assert graph is not accepted
    assert report["selected_source"] == "beta"
    assert [row["source"] for row in report["candidates"]] == ["beta", "zeta", "alpha"]
    candidate = DepictionCandidate("worker", accepted, 1.5, metadata={"origin": "worker"})
    assert candidate.source == "worker"
    assert candidate.graph is accepted
    assert candidate.score == 1.5
    assert not candidate.rejected
    assert candidate.rejection_reason == ""
    assert candidate.metadata == {"origin": "worker"}


def test_empty_smiles_and_worker_error_keep_existing_errors(monkeypatch) -> None:
    with pytest.raises(ValueError, match="SMILES vacío"):
        rdkit_io.smiles_to_depiction_candidates("  ")

    monkeypatch.setattr(rdkit_safe, "smiles_depict_candidates_isolated", lambda *_args, **_kwargs: ([], "worker_down"))
    with pytest.raises(RuntimeError, match="RDKit worker no disponible"):
        rdkit_io.smiles_to_depiction_candidates("CC")


def test_molblock_parse_rejection_preserves_source_and_metadata(monkeypatch) -> None:
    monkeypatch.setattr(
        rdkit_safe,
        "smiles_depict_candidates_isolated",
        lambda *_args, **_kwargs: ([{"source": "broken", "ok": True, "molblock": "bad", "metadata": {"raw": 1}}], None),
    )
    monkeypatch.setattr(rdkit_io, "molfile_to_molgraph", lambda _molblock: (_ for _ in ()).throw(ValueError("invalid_molblock")))

    candidates = rdkit_io.smiles_to_depiction_candidates("CC")

    assert len(candidates) == 1
    assert candidates[0].source == "broken"
    assert candidates[0].rejected
    assert candidates[0].rejection_reason == "invalid_molblock"
    assert candidates[0].metadata == {"raw": 1}


def test_vancomycin_smiles_import_produces_non_degenerate_graph() -> None:
    try:
        graph = rdkit_io.smiles_to_molgraph_best_depiction(VANCOMYCIN_SMILES, timeout_s=20.0)
    except RuntimeError as exc:
        pytest.skip(f"Worker RDKit no disponible en este entorno: {exc}")

    assert len(graph.atoms) > 80
    assert len(graph.bonds) > 80
    lengths = [
        math.hypot(graph.atoms[bond.a2_id].x - graph.atoms[bond.a1_id].x, graph.atoms[bond.a2_id].y - graph.atoms[bond.a1_id].y)
        for bond in graph.bonds.values()
    ]
    assert sum(1 for length in lengths if length > 10.0) / len(lengths) >= 0.70
    xs = [atom.x for atom in graph.atoms.values()]
    ys = [atom.y for atom in graph.atoms.values()]
    assert max(xs) - min(xs) > 200.0
    assert max(ys) - min(ys) > 200.0
    assert any(abs(atom.x) > 1e-3 or abs(atom.y) > 1e-3 for atom in graph.atoms.values())


def test_vancomycin_depiction_prefers_coordgen_when_available() -> None:
    try:
        candidates = rdkit_io.smiles_to_depiction_candidates(VANCOMYCIN_SMILES, timeout_s=20.0)
    except RuntimeError as exc:
        pytest.skip(f"Worker RDKit no disponible en este entorno: {exc}")
    coordgen = [candidate for candidate in candidates if candidate.source == "rdcoordgen" and not candidate.rejected]
    if not coordgen:
        pytest.skip("RDKit CoordGen no disponible en este entorno")
    compute2d = [candidate for candidate in candidates if candidate.source == "rddepictor_compute2d" and not candidate.rejected]
    if not compute2d:
        pytest.skip("Compute2D no produjo candidato comparable")
    ordered_sources = [candidate.source for candidate in candidates if not candidate.rejected]
    if coordgen[0].score <= compute2d[0].score:
        assert ordered_sources.index("rdcoordgen") < ordered_sources.index("rddepictor_compute2d")


def test_compute2d_worker_still_handles_simple_smiles() -> None:
    try:
        graph = rdkit_io.smiles_to_molgraph_best_depiction("CCO", timeout_s=10.0)
    except RuntimeError as exc:
        pytest.skip(f"Worker RDKit no disponible en este entorno: {exc}")
    assert len(graph.atoms) >= 3
    lengths = [
        math.hypot(graph.atoms[bond.a2_id].x - graph.atoms[bond.a1_id].x, graph.atoms[bond.a2_id].y - graph.atoms[bond.a1_id].y)
        for bond in graph.bonds.values()
    ]
    assert lengths and min(lengths) > 1e-3


def test_depiction_score_penalizes_circular_donut_layout_for_block_rich_graph() -> None:
    circular = _block_rich_graph(circular=True)
    elongated = _block_rich_graph(circular=False)

    circular_score, _ = score_imported_depiction(circular, target_bond_length=40.0)
    elongated_score, _ = score_imported_depiction(elongated, target_bond_length=40.0)

    assert circular_score > elongated_score


def _block_rich_graph(*, circular: bool) -> MolGraph:
    graph = MolGraph()
    rings: list[list[int]] = []
    count = 5
    for idx in range(count):
        if circular:
            angle = 2.0 * math.pi * idx / count
            cx, cy = 120.0 * math.cos(angle), 120.0 * math.sin(angle)
        else:
            cx, cy = idx * 130.0, 35.0 * math.sin(idx)
        rings.append(_add_hexagon(graph, cx, cy))
    for idx in range(count - 1):
        _connect(graph, rings[idx][0], rings[idx + 1][3])
    return graph


def _add_hexagon(graph: MolGraph, cx: float, cy: float) -> list[int]:
    atoms = []
    for idx in range(6):
        angle = math.radians(idx * 60.0)
        atoms.append(graph.add_atom("C", cx + 40.0 * math.cos(angle), cy + 40.0 * math.sin(angle)).id)
    for idx in range(6):
        graph.add_bond(atoms[idx], atoms[(idx + 1) % 6], order=1, is_aromatic=True)
    return atoms


def _connect(graph: MolGraph, left: int, right: int) -> None:
    lx, ly = graph.atoms[left].x, graph.atoms[left].y
    rx, ry = graph.atoms[right].x, graph.atoms[right].y
    mid = graph.add_atom("C", (lx + rx) * 0.5, (ly + ry) * 0.5).id
    graph.add_bond(left, mid, order=1)
    graph.add_bond(mid, right, order=1)
