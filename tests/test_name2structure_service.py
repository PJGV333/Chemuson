"""Pruebas de Name→Structure desacoplado."""

from __future__ import annotations

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.name2structure import (
    NameToStructureResult,
    StaticNameConnector,
    resolve_name_to_structure,
)
import chemuson.name2structure.service as service


def test_static_connector_resolves_common_name_without_network(monkeypatch) -> None:
    def fake_smiles_to_graph(smiles: str, timeout_s: float):
        graph = MolGraph()
        graph.add_atom("O" if smiles == "O" else "C", 0.0, 0.0)
        return graph, ""

    monkeypatch.setattr(service, "_smiles_to_graph", fake_smiles_to_graph)

    result = resolve_name_to_structure("agua", allow_network=False)

    assert result.ok
    assert result.source == "offline-common"
    assert result.smiles == "O"
    assert result.confidence > 0.0
    assert result.graph is not None
    assert len(result.graph.atoms) == 1


def test_resolver_uses_connectors_in_order() -> None:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)

    class MissingConnector:
        source = "missing"

        def resolve(self, name: str, timeout_s: float = 8.0):
            return NameToStructureResult(name, None, self.source, 0.0, message="not_found")

    class HitConnector:
        source = "hit"

        def resolve(self, name: str, timeout_s: float = 8.0):
            return NameToStructureResult(name, graph, self.source, 0.95, smiles="C")

    result = resolve_name_to_structure(
        "methane",
        allow_network=False,
        connectors=[MissingConnector(), HitConnector()],
    )

    assert result.ok
    assert result.source == "hit"
    assert result.confidence == 0.95


def test_static_connector_reports_not_found_without_rdkit_call(monkeypatch) -> None:
    def fail_smiles_to_graph(smiles: str, timeout_s: float):
        raise AssertionError("should not be called")

    monkeypatch.setattr(service, "_smiles_to_graph", fail_smiles_to_graph)

    result = StaticNameConnector(entries={}).resolve("unknown")

    assert not result.ok
    assert result.message == "not_found"


def test_static_connector_uses_internal_fallback_when_rdkit_worker_fails(monkeypatch) -> None:
    import chemuson.chemio.rdkit_safe as rdkit_safe

    def fail_worker(_smiles: str, timeout_s: float):
        return None, "rdkit_unavailable"

    monkeypatch.setattr(rdkit_safe, "smiles_to_molgraph_isolated", fail_worker)

    result = resolve_name_to_structure("ethanol", allow_network=False)

    assert result.ok
    assert result.source == "offline-common"
    assert result.smiles == "CCO"
    assert result.graph is not None
    assert [atom.element for atom in result.graph.atoms.values()] == ["C", "C", "O"]
