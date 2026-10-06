from __future__ import annotations

from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated
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
        resolver=lambda name: _reference("CCO", name),
    )

    assert result.status is MolecularIdentityStatus.VERIFIED
    assert result.reference_identifier == "offline-fake-reference:ethanol"


def test_valid_but_different_structure_is_a_mismatch():
    result = verify_molecular_identity(
        "Dibuja la cafeína",
        _graph("CCN"),
        resolver=lambda name: _reference("CCO", name),
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
        "Draw caffeine", _graph("CCO"), resolver=lambda _name: missing
    ).status is MolecularIdentityStatus.UNVERIFIED
    assert verify_molecular_identity(
        "Draw caffeine", _graph("CCO"), resolver=lambda _name: unavailable
    ).status is MolecularIdentityStatus.REFERENCE_ERROR
    assert verify_molecular_identity(
        "Draw caffeine",
        _graph("CCO"),
        resolver=lambda _name: (_ for _ in ()).throw(RuntimeError("private")),
    ).reason_code == "resolver_error"


def test_open_ended_prompt_does_not_call_reference_resolver():
    calls = []
    result = verify_molecular_identity(
        "Generate a molecule with three fused rings",
        _graph("CCO"),
        resolver=lambda name: calls.append(name),
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
        resolver=lambda name: calls.append(name)
        or _reference(CHOLESTEROL_REFERENCE_SMILES, name),
        allow_network=False,
    )

    assert calls == ["colesterol"]
    assert len(cholesterol.atoms) == 28
    assert result.status is MolecularIdentityStatus.MISMATCH
    assert result.requested_name == "colesterol"


def test_identity_result_depends_on_requested_name_and_graph_not_model_narrative():
    reference = _reference("CCO", "ethanol")
    graph = _graph("CCN")

    def resolver(_name: str) -> NameToStructureResult:
        return reference

    first = verify_molecular_identity("Draw ethanol", graph, resolver=resolver)
    second = verify_molecular_identity("Draw ethanol", graph, resolver=resolver)

    assert first == second
    assert first.status is MolecularIdentityStatus.MISMATCH
