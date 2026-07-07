from __future__ import annotations

from chemuson.clean2d import (
    Clean2DCandidate,
    Clean2DMode,
    Clean2DResult,
    stable_clean2d_rejection_reason,
)
from chemuson.gui.controllers.clean2d_controller import Clean2DAttempt


def test_stable_rejection_reason_mapping_preserves_contract_vocabulary() -> None:
    cases = {
        "coordenada_no_finita:7": "invalid-coordinates",
        "cambio_metadatos_atomos": "invariant-violation",
        "cambia_estereoquimica": "stereo-risk",
        "reparacion_rechazada_por_integridad_de_enlaces": "boundary-bond-risk",
        "nuevos_cruces_enlaces:2": "new-crossing-risk",
        "colision_no_enlazada": "collision-risk",
        "anillo_colapsado": "collapsed-ring-risk",
        "desplazamiento_maximo_excesivo:999": "excessive-displacement",
        "candidato_no_mejora_suficiente": "worse-quality",
        "rdkit_aislado_no_disponible:test": "backend-failure",
    }

    for internal_reason, stable_reason in cases.items():
        assert stable_clean2d_rejection_reason(internal_reason) == stable_reason


def test_clean2d_result_exposes_stable_operation_states() -> None:
    before = {1: (0.0, 0.0), 2: (40.0, 0.0)}
    current = Clean2DCandidate(source="current", coords=before, novelty=0.0)
    applied = Clean2DCandidate(source="clean2d_v2", coords={1: (0.0, 1.0), 2: (40.0, 1.0)}, novelty=1.0)
    preserve = Clean2DCandidate(
        source="current",
        coords=before,
        metadata={"complex_preserve_only": True},
    )
    rejected = Clean2DCandidate(
        source="rdkit_isolated",
        coords={},
        rejected=True,
        rejection_reason="rdkit_aislado_no_disponible:test",
    )

    assert _result(current, before).result_state == "no-op"
    assert _result(applied, before).result_state == "applied"
    assert _result(preserve, before).result_state == "preserve-only"

    failed = Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids={1, 2},
        before_coords=before,
        candidates=(),
        selected=None,
        rejected=(rejected,),
        reason=rejected.rejection_reason,
    )
    assert failed.result_state == "failed-controlled"
    assert failed.stable_reason == "backend-failure"
    assert failed.candidate_sources == ("rdkit_isolated",)


def test_clean2d_candidate_keeps_backend_detail_and_stable_reason() -> None:
    candidate = Clean2DCandidate(
        source="rdkit_isolated",
        coords={},
        rejected=True,
        rejection_reason="rdkit_aislado_no_disponible:worker_timeout",
    )

    assert candidate.outcome_state == "rejected"
    assert candidate.rejection_reason.endswith("worker_timeout")
    assert candidate.stable_rejection_reason == "backend-failure"


def test_controller_attempt_exposes_stable_diagnostic_reason() -> None:
    attempt = Clean2DAttempt(
        rejected=True,
        message="Limpieza 2D omitida: reparación rechazada por integridad",
        details={"stable_rejection_reason": "stereo-risk", "backend_detail": "signature mismatch"},
    )

    assert attempt.result_state == "rejected"
    assert attempt.stable_reason == "stereo-risk"
    assert attempt.details["backend_detail"] == "signature mismatch"


def _result(candidate: Clean2DCandidate, before: dict[int, tuple[float, float]]) -> Clean2DResult:
    return Clean2DResult(
        mode=Clean2DMode.QUICK,
        atom_ids=set(before),
        before_coords=before,
        candidates=(candidate,),
        selected=candidate,
    )
