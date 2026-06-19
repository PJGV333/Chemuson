from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.chemio import rdkit_io, rdkit_safe
from chemuson.chemio.rdkit_safe import rdkit_worker_diagnostics, smiles_to_molgraph_isolated
from chemuson.core.model import MolGraph


def _simple_graph() -> MolGraph:
    graph = MolGraph()
    a1 = graph.add_atom("C", 0.0, 0.0)
    a2 = graph.add_atom("C", 40.0, 0.0)
    graph.add_bond(a1.id, a2.id, order=1)
    return graph


def test_smiles_to_molgraph_prefers_isolated_worker_when_direct_rdkit_disabled(monkeypatch) -> None:
    monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: False)
    monkeypatch.setattr(rdkit_safe, "smiles_to_molgraph_isolated", lambda _smiles, timeout_s=8.0: (_simple_graph(), None))

    graph = rdkit_io.smiles_to_molgraph("CC")

    assert len(graph.atoms) >= 2


def test_smiles_to_molgraph_reports_worker_diagnostics_on_failure(monkeypatch) -> None:
    monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: False)
    monkeypatch.setattr(
        rdkit_safe,
        "smiles_to_molgraph_isolated",
        lambda _smiles, timeout_s=8.0: (None, "rdkit_unavailable: test"),
    )

    with pytest.raises(RuntimeError) as excinfo:
        rdkit_io.smiles_to_molgraph("CC")

    message = str(excinfo.value)
    assert "RDKit worker no disponible" in message
    assert "rdkit_unavailable" in message
    assert "sys.executable" in message or sys.executable in message


def test_isolated_smiles_import_generates_non_degenerate_coordinates() -> None:
    graph, error = smiles_to_molgraph_isolated("CC", timeout_s=5.0)
    if graph is None:
        pytest.skip(f"Worker RDKit no disponible en este entorno: {error}")

    assert len(graph.atoms) >= 2
    lengths = [
        math.hypot(
            graph.atoms[bond.a2_id].x - graph.atoms[bond.a1_id].x,
            graph.atoms[bond.a2_id].y - graph.atoms[bond.a1_id].y,
        )
        for bond in graph.bonds.values()
    ]
    assert lengths
    assert max(lengths) > 1e-3


def test_rdkit_worker_diagnostics_returns_environment_payload() -> None:
    payload = rdkit_worker_diagnostics(timeout_s=5.0)

    assert isinstance(payload, dict)
    assert "python_executable" in payload
    if payload.get("ok"):
        assert "rdkit_version" in payload
