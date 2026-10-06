from __future__ import annotations

from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
import chemuson.name2structure.identity as identity_module
from chemuson.name2structure import (
    MolecularIdentityStatus,
    NameToStructureResult,
    extract_requested_molecule_name,
    verify_molecular_identity,
)


CHOLESTEROL_REFERENCE_SMILES = (
    "CC(C)CCC[C@H](C)[C@@H]1CC[C@@H]2[C@@H]3CC=C4"
    "C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
)


def _graph(smiles: str):
    graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=5.0)
    assert error is None
    assert graph is not None
    return graph


def _reference(smiles: str, name: str) -> NameToStructureResult:
    return NameToStructureResult(
        query=name,
        graph=_graph(smiles),
        source="offline-fake-reference",
        confidence=1.0,
        smiles=smiles,
        resolved_name=name,
    )


def test_equivalent_structures_with_different_smiles_are_verified():
    result = verify_molecular_identity(
        "Draw ethanol",
        _graph("OCC"),
        resolver=lambda name, *, allow_network: _reference("CCO", name),
    )

    assert result.status is MolecularIdentityStatus.VERIFIED
    assert result.reference_identifier == "offline-fake-reference:ethanol"


def test_valid_but_different_structure_is_a_mismatch():
    result = verify_molecular_identity(
        "Dibuja la cafeína",
        _graph("CCN"),
        resolver=lambda name, *, allow_network: _reference("CCO", name),
    )

    assert result.status is MolecularIdentityStatus.MISMATCH
    assert result.requested_name == "cafeína"


def test_missing_and_failing_reference_are_not_verified():
    missing = NameToStructureResult(
        query="caffeine",
        graph=None,
        source="fake",
        confidence=0.0,
        message="not_found",
    )
    unavailable = NameToStructureResult(
        query="caffeine",
        graph=None,
        source="fake",
        confidence=0.0,
        message="TimeoutError",
    )

    assert verify_molecular_identity(
        "Draw caffeine",
        _graph("CCO"),
        resolver=lambda _name, *, allow_network: missing,
    ).status is MolecularIdentityStatus.UNVERIFIED
    assert verify_molecular_identity(
        "Draw caffeine",
        _graph("CCO"),
        resolver=lambda _name, *, allow_network: unavailable,
    ).status is MolecularIdentityStatus.REFERENCE_ERROR
    assert verify_molecular_identity(
        "Draw caffeine",
        _graph("CCO"),
        resolver=lambda _name, *, allow_network: (_ for _ in ()).throw(RuntimeError("private")),
    ).reason_code == "resolver_error"


def test_open_ended_prompt_does_not_call_reference_resolver():
    calls = []
    result = verify_molecular_identity(
        "Generate a molecule with three fused rings",
        _graph("CCO"),
        resolver=lambda name, *, allow_network: calls.append(name),
    )

    assert result.status is MolecularIdentityStatus.NOT_APPLICABLE
    assert calls == []
    assert extract_requested_molecule_name("Dibuja el colesterol") == "colesterol"
    assert extract_requested_molecule_name("Draw caffeine") == "caffeine"
    assert extract_requested_molecule_name("Genera una molécula con tres anillos") is None


def test_cholesterol_semantic_mismatch_01_uses_offline_reference():
    cholesterol = _graph(CHOLESTEROL_REFERENCE_SMILES)
    proposal = _graph("CCO")
    calls = []

    result = verify_molecular_identity(
        "Dibuja colesterol",
        proposal,
        resolver=lambda name, *, allow_network: calls.append((name, allow_network))
        or _reference(CHOLESTEROL_REFERENCE_SMILES, name),
        allow_network=False,
    )

    assert calls == [("colesterol", False)]
    assert len(cholesterol.atoms) == 28
    assert result.status is MolecularIdentityStatus.MISMATCH
    assert result.requested_name == "colesterol"


def test_identity_defaults_to_offline_and_a_missing_reference_stays_unverified(monkeypatch):
    missing = NameToStructureResult(
        query="unknown molecule",
        graph=None,
        source="offline-common",
        confidence=0.0,
        message="not_found",
    )
    calls = []

    def resolve(name, *, allow_network, timeout_s):
        calls.append((name, allow_network, timeout_s))
        return missing

    monkeypatch.setattr(identity_module, "resolve_name_to_structure", resolve)
    result = verify_molecular_identity("Draw unknown molecule", _graph("CCO"))

    assert calls == [("unknown molecule", False, 8.0)]
    assert result.status is MolecularIdentityStatus.UNVERIFIED
    assert result.reason_code == "reference_not_found_offline"


def test_offline_reference_can_verify_identity_without_external_lookup():
    calls = []

    def resolve(name, *, allow_network):
        calls.append((name, allow_network))
        return _reference("CCO", name)

    result = verify_molecular_identity("Draw ethanol", _graph("OCC"), resolver=resolve)

    assert calls == [("ethanol", False)]
    assert result.status is MolecularIdentityStatus.VERIFIED


def test_external_reference_requires_explicit_opt_in_and_propagates_policy():
    calls = []

    def resolve(name, *, allow_network):
        calls.append((name, allow_network))
        return _reference("CCO", name)

    result = verify_molecular_identity(
        "Draw ethanol",
        _graph("CCO"),
        resolver=resolve,
        allow_network=True,
    )

    assert calls == [("ethanol", True)]
    assert result.status is MolecularIdentityStatus.VERIFIED


def test_disabled_identity_verification_skips_resolution_and_canonicalization(monkeypatch):
    def unexpected(*_args, **_kwargs):
        raise AssertionError("disabled verification must do no work")

    monkeypatch.setattr(identity_module, "resolve_name_to_structure", unexpected)
    from chemuson.chemio import rdkit_safe

    monkeypatch.setattr(rdkit_safe, "molgraph_to_inchi_isolated", unexpected)
    result = verify_molecular_identity(
        "Draw ethanol",
        _graph("CCO"),
        enabled=False,
        allow_network=True,
    )

    assert result.status is MolecularIdentityStatus.NOT_APPLICABLE
    assert result.reason_code == "verification_disabled"


def test_identity_result_depends_on_requested_name_and_graph_not_model_narrative():
    reference = _reference("CCO", "ethanol")
    graph = _graph("CCN")

    def resolver(
        _name: str, *, allow_network: bool
    ) -> NameToStructureResult:
        assert allow_network is False
        return reference

    first = verify_molecular_identity("Draw ethanol", graph, resolver=resolver)
    second = verify_molecular_identity("Draw ethanol", graph, resolver=resolver)

    assert first == second
    assert first.status is MolecularIdentityStatus.MISMATCH
