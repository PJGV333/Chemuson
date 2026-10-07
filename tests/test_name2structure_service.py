"""Pruebas de Name→Structure desacoplado."""

from __future__ import annotations

import json

import pytest

from chemuson.core.model import MolGraph
from chemuson.name2structure import (
    NameToStructureResult,
    PubChemNameConnector,
    StaticNameConnector,
    extract_requested_molecule_name,
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


def test_offline_resolver_uses_pubchem_cache_without_external_fetch(monkeypatch, tmp_path) -> None:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    cache_path = tmp_path / ".chemuson" / "name2structure_cache.json"
    cache_path.parent.mkdir()
    cache_path.write_text(
        json.dumps(
            {
                "cached molecule": {
                    "smiles": "C",
                    "resolved_name": "cached molecule",
                    "confidence": 0.86,
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(service.Path, "home", lambda: tmp_path)
    monkeypatch.setattr(service, "_smiles_to_graph", lambda *_args, **_kwargs: (graph, ""))
    fetches = []
    monkeypatch.setattr(
        PubChemNameConnector,
        "_fetch_smiles",
        lambda *_args, **_kwargs: fetches.append(True),
    )

    result = resolve_name_to_structure("cached molecule", allow_network=False)

    assert result.ok
    assert result.source == "pubchem"
    assert result.from_cache is True
    assert result.graph is graph
    missing = resolve_name_to_structure("not in cache", allow_network=False)
    assert not missing.ok
    assert missing.message == "not_found"
    assert fetches == []


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


class _FakePugResponse:
    def __init__(self, payload: dict[str, object]) -> None:
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *_args) -> None:
        return None

    def read(self) -> bytes:
        return json.dumps(self._payload).encode("utf-8")


def _mock_pug_response(monkeypatch, payload: dict[str, object], calls=None) -> None:
    def fake_urlopen(request, timeout):
        if calls is not None:
            calls.append((request.full_url, timeout))
        return _FakePugResponse(payload)

    monkeypatch.setattr(service, "urlopen", fake_urlopen)


def test_pubchem_current_pug_contract_requests_and_prefers_smiles(monkeypatch, tmp_path) -> None:
    payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 1,
                    "SMILES": "C[C@H](O)N",
                    "ConnectivitySMILES": "CC(O)N",
                    "IUPACName": "2-aminopropan-1-ol",
                }
            ]
        }
    }
    calls = []
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    _mock_pug_response(monkeypatch, payload, calls)
    monkeypatch.setattr(service, "_smiles_to_graph", lambda smiles, timeout_s: (graph, ""))

    result = PubChemNameConnector(cache_path=tmp_path / "cache.json").resolve("example")

    assert result.ok
    assert result.smiles == "C[C@H](O)N"
    assert result.resolved_name == "2-aminopropan-1-ol"
    assert len(calls) == 1
    requested_url, timeout = calls[0]
    assert "/property/SMILES,ConnectivitySMILES,IUPACName/JSON" in requested_url
    assert "IsomericSMILES,CanonicalSMILES" not in requested_url
    assert timeout >= 1.0


def test_pubchem_uses_connectivity_smiles_when_smiles_is_absent(monkeypatch, tmp_path) -> None:
    payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 2,
                    "ConnectivitySMILES": "CCO",
                    "IUPACName": "ethanol",
                }
            ]
        }
    }
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    _mock_pug_response(monkeypatch, payload)
    monkeypatch.setattr(service, "_smiles_to_graph", lambda smiles, timeout_s: (graph, ""))

    result = PubChemNameConnector(cache_path=tmp_path / "cache.json").resolve("ethanol")

    assert result.ok
    assert result.smiles == "CCO"


@pytest.mark.parametrize(
    ("legacy_property", "legacy_smiles"),
    [("IsomericSMILES", "C[C@H](O)N"), ("CanonicalSMILES", "CC(O)N")],
)
def test_pubchem_parser_retains_legacy_smiles_compatibility(
    monkeypatch, tmp_path, legacy_property, legacy_smiles
) -> None:
    payload = {
        "PropertyTable": {
            "Properties": [
                {"CID": 3, legacy_property: legacy_smiles, "IUPACName": "example"}
            ]
        }
    }
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    _mock_pug_response(monkeypatch, payload)
    monkeypatch.setattr(service, "_smiles_to_graph", lambda smiles, timeout_s: (graph, ""))

    result = PubChemNameConnector(cache_path=tmp_path / "cache.json").resolve("legacy")

    assert result.ok
    assert result.smiles == legacy_smiles


def test_pubchem_missing_smiles_returns_controlled_empty_smiles(monkeypatch, tmp_path) -> None:
    _mock_pug_response(
        monkeypatch,
        {"PropertyTable": {"Properties": [{"CID": 4, "IUPACName": "example"}]}},
    )

    result = PubChemNameConnector(cache_path=tmp_path / "cache.json").resolve("empty")

    assert not result.ok
    assert result.graph is None
    assert result.message == "empty_smiles"


def test_pubchem_smiles_rejected_by_chemio_is_never_a_reference(monkeypatch, tmp_path) -> None:
    cache_path = tmp_path / "cache.json"
    _mock_pug_response(
        monkeypatch,
        {
            "PropertyTable": {
                "Properties": [
                    {"CID": 5, "SMILES": "C1CC", "ConnectivitySMILES": "C1CC"}
                ]
            }
        },
    )

    result = PubChemNameConnector(cache_path=cache_path).resolve("invalid")

    assert not result.ok
    assert result.graph is None
    assert result.smiles == "C1CC"
    assert result.message
    assert not cache_path.exists()


@pytest.mark.parametrize(
    ("prompt", "original_query", "canonical_query"),
    [
        ("Dibuja la tetrandrina", "tetrandrina", "tetrandrine"),
        ("Dibuja colesterol", "colesterol", "cholesterol"),
    ],
)
def test_verified_name_aliases_preserve_original_and_canonical_queries(
    prompt, original_query, canonical_query
) -> None:
    requested_name = extract_requested_molecule_name(prompt)
    assert requested_name == original_query
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    calls = []

    class RecordingConnector:
        source = "pubchem"

        def resolve(self, name: str, timeout_s: float = 8.0):
            calls.append(name)
            return NameToStructureResult(
                name,
                graph,
                self.source,
                0.86,
                smiles="CN1CC",
                resolved_name=f"{canonical_query} IUPAC",
            )

    result = resolve_name_to_structure(
        requested_name,
        allow_network=False,
        connectors=[RecordingConnector()],
    )

    assert calls == [canonical_query]
    assert result.query == original_query
    assert result.resolved_query == canonical_query
    assert result.resolved_name == f"{canonical_query} IUPAC"
    assert result.source == "pubchem"
