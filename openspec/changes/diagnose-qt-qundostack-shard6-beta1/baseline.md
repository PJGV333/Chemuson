# Baseline — SIGABRT Qt/QUndoStack shard 6

Fecha: 2026-10-21 UTC. Rama `stabilize/v0.3.0-beta.1`, HEAD limpio `aff565a064f17dccf98b0b4f23c7ddeba593ca31`, remoto igual. Base de comparación pedida: `61f7647b683db3eca19ea06e5ec82ab2897d640a`. No se ha modificado código desde el commit/push.

## Baseline del SHA actual

- `chem/bin/python`: Python 3.14.8, PyQt6 6.11.0, Qt runtime 6.11.2, RDKit 2026.03.6, pytest 9.1.1.
- `python -m compileall src tests tools packaging`: exit 0.
- `QT_QPA_PLATFORM=offscreen chem/bin/python -m pytest --collect-only -q`: **2.598 collected**, exit 0.
- `timeout --signal=TERM --kill-after=10s 480s chem/bin/python -m pytest -q`: exit 124 a 480 s (~30%); no summary. No se trata como PASS ni como reproducción del SIGABRT. Log temporal exacto: `baseline-current/pytest-full.txt` en el área de diagnóstico local.
- Ruff `F401,F811,F821,E722,E741`: exit 1 únicamente por F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- El planner actual asigna 2.598 IDs, distribución `[248,243,327,362,449,268,373,328]`, SHA manifest `b2cb90fdc74a1c2417885557f9391f12924235c3043d2e1a73fee0340027248f`.

## Evidencia del fallo actual

- Actions run [#39689871971](https://github.com/PJGV333/Chemuson/actions/runs/39689871971), HEAD exacto `aff565a064f17dccf98b0b4f23c7ddeba593ca31`: `pytest-plan` PASS, shards 0–5 y 7 PASS, shard 6 **SIGABRT** sin JUnit; el test iniciado es `tests/test_ui_polish.py::test_onboarding_close_releases_overlay`. Smoke `windows-smoke` PASS e `flatpak-smoke` PASS. El verifier y workflow concluyen failure. Artefactos fueron consultados desde el run; job shard 6 `119533389359`, Windows `119533119387`, Flatpak `119533119538`.
- Evidencia local previa a este OpenSpec: shard 6 de 373 IDs abortó dos veces (144.659 s y 143.567 s), en la misma línea/test con `QUndoStack::~QUndoStack`/`cleanChanged` en el stack nativo; 308 timings completos, JUnit ausente, no timeout. `tests/test_ui_polish.py` aislado: 38 passed. JUnit/timings/plan actuales se conservan temporalmente en `/tmp/chemuson-final-stabilization-shards/`.
- Actions base [#39679607676](https://github.com/PJGV333/Chemuson/actions/runs/39679607676), SHA `61f7647b683db3eca19ea06e5ec82ab2897d640a`: 2.597 IDs y ejecución normal de 306/306 en shard 6 (304 passed, 2 failed, 0 errors, sin señal fatal). La misma prueba `test_onboarding_close_releases_overlay` sí estaba asignada a ese shard y pasó. Windows y Flatpak smoke también pasaron. El conjunto anterior de IDs/orden del shard 6 no era igual al actual; esto no basta aún para excluir regresión.
## Comparación causal A/B (sin cambios)

- Se extrajeron snapshots `git archive` de `61f7647b683db3eca19ea06e5ec82ab2897d640a` y `aff565a064f17dccf98b0b4f23c7ddeba593ca31` dentro del área temporal del repo; no se crearon worktrees. En ambos procesos `chemuson.__file__` resolvió al snapshot correspondiente.
- Runtime idéntico: Python 3.11.15, pytest 9.1.1, PyQt6 6.11.0 / Qt 6.11.2, RDKit 2026.09.1; `QT_QPA_PLATFORM=offscreen`. Se reprodujeron los mismos **372 IDs en el mismo orden** del plan actual, excluyendo sólo el nuevo node `test_compchem_worker_stops_its_thread_without_main_event_loop` que no existe en la base.
- Sin parche, base y candidato abortan con exit 134 en `tests/test_ui_polish.py::test_onboarding_close_releases_overlay`, ambos tras 307 resultados temporizados (**305 passed, 2 skipped**), sin JUnit. La ruta CompChem anterior pasa en ambos. El stack fatal incluye `QUndoStack::~QUndoStack`/`cleanChanged` durante `sendPostedEvents`.
- Experimento red/green temporal, cambiando sólo el receiver de `QApplication.sendPostedEvents(None, DeferredDelete)` a `QApplication.sendPostedEvents(overlay, DeferredDelete)`: base exit 0, **370 passed, 2 skipped**, 372 timings, JUnit completo, 205.16 s; candidato exit 0 con los mismos resultados, 206.01 s. Las fuentes temporales de estos experimentos se eliminarán; logs/JUnit/timings se preservan como evidencia `evidence-*-py311-*` junto a este documento.
- Diagnóstico: el test pedía destruir un overlay pero vaciaba todos los DeferredDelete globales, incluyendo widgets ajenos acumulados por el orden del shard. Es un defecto de aislamiento del test/lifecycle ordenado, preexistente en ambos SHAs; no hay evidencia que justifique cambiar `CanvasTabManager` ni otra producción.

El arreglo autorizado se limita al receiver de `sendPostedEvents` en `tests/test_ui_polish.py`, más un comentario que explique el scope. No se alteran el planner, asserts, skips ni timeouts.
