# Validación — SIGABRT Qt/QUndoStack del shard 6

Fecha: 2026-10-21 UTC. Rama `stabilize/v0.3.0-beta.1`. Esta validación documenta el diagnóstico y el fix local; los gates remotos pendientes se anotarán al validar el SHA final. No incluye auditoría de paquetes ni aceptación manual de la beta.

## Causa y comparación controlada

El fallo no era una regresión de `aff565a` frente a `61f7647b`. El test `test_onboarding_close_releases_overlay` solicitaba la destrucción de un overlay, pero `QApplication.sendPostedEvents(None, DeferredDelete)` descargaba todos los DeferredDelete de la aplicación, incluidos widgets/undo stacks ajenos acumulados por el orden del shard. La destrucción diferida ajena podía abortar en `QUndoStack::~QUndoStack`/`cleanChanged`. La corrección restringe el envío al receiver `overlay`; no modifica producción ni cambia la aserción `sip.isdeleted(overlay)`.

Evidencia A/B en snapshots temporales de `git archive`, sin worktrees, con los mismos 372 IDs del shard 6 en el orden actual y excluyendo sólo el node nuevo de CompChem ausente en base:

| Código fuente | Test sin corregir | Sólo receiver del overlay corregido |
|---|---:|---:|
| Base `61f7647b683db3eca19ea06e5ec82ab2897d640a` | exit 134, SIGABRT; 307 timings completos (305 passed, 2 skipped), JUnit ausente | exit 0; 370 passed, 2 skipped; 372 timings y JUnit completos; 205.16 s |
| Candidato `aff565a064f17dccf98b0b4f23c7ddeba593ca31` | exit 134, SIGABRT en el mismo test/línea; 307 timings completos (305 passed, 2 skipped), JUnit ausente | exit 0; 370 passed, 2 skipped; 372 timings y JUnit completos; 206.01 s |

Runtime A/B idéntico: Python 3.11.15, pytest 9.1.1, PyQt6 6.11.0 / Qt 6.11.2, RDKit 2026.09.1, Linux y `QT_QPA_PLATFORM=offscreen`; se verificó `chemuson.__file__` en cada snapshot. Los logs, JUnit, timings y exits A/B están preservados en `evidence-*`; sus hashes están en `evidence-sha256.txt`. La prueba `test_ui_polish.py` aislada también pasa (38 tests). La secuencia base/candidato pasa el test CompChem anterior en ambos; el aborto aparece en onboarding bajo el orden del shard.

## Cambio

- `tests/test_ui_polish.py`: `sendPostedEvents` recibe el overlay concreto. Un comentario explica por qué la prueba no debe procesar deletes pendientes de otros tests.
- No se tocaron `src/`, callbacks `CanvasTabManager`, timeouts, skips, planner, manifiesto ni `KNOWN_BASELINE_EXCEPTIONS`.
- La colección permanece en **2.598 IDs**; no apareció node nuevo ni fue necesario modificar `tests/ci/expected_pytest_nodeids.txt`.

## Verificación local del checkout

Runtime principal: Python 3.14.8, pytest 9.1.1, Qt 6.11.2, RDKit 2026.03.6, `QT_QPA_PLATFORM=offscreen`.

- Focal UI + contrato del workflow: **42 passed** en 16.95 s.
- `tests/architecture`: **283 passed** en 9.85 s.
- `python -m compileall -q src tests tools packaging`: PASS.
- Ruff focal `tests/test_ui_polish.py --select F401,F811,F821,E722,E741`: PASS.
- Ruff global con el mismo selector: únicamente el F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se cambió.
- `openspec validate diagnose-qt-qundostack-shard6-beta1 --strict`, `git diff --check`: PASS.
- Plan CI: 2.598 IDs; asignaciones `[248, 243, 327, 362, 449, 268, 373, 328]`; hash esperado del set `b2cb90fdc74a1c2417885557f9391f12924235c3043d2e1a73fee0340027248f`.
- Ocho shards ejecutados **secuencialmente**, cada uno con cap de 300 s; los ocho exits son 0. JUnit y JSONL de timings cubren los IDs asignados exactamente. El verificador local reporta:

```text
shard 0: exit=0, tests=248/248, timed=248/248, elapsed=87.662s
shard 1: exit=0, tests=243/243, timed=243/243, elapsed=77.634s
shard 2: exit=0, tests=327/327, timed=327/327, elapsed=68.128s
shard 3: exit=0, tests=362/362, timed=362/362, elapsed=86.170s
shard 4: exit=0, tests=449/449, timed=449/449, elapsed=60.666s
shard 5: exit=0, tests=268/268, timed=268/268, elapsed=179.352s
shard 6: exit=0, tests=373/373, timed=373/373, elapsed=207.851s
shard 7: exit=0, tests=328/328, timed=328/328, elapsed=68.897s
pytest shard coverage verified: 2598 unique node IDs, 8 complete successful shards
```

Totales JUnit: **2.576 passed, 22 skipped, 0 failures, 0 errors**. `verify.log` y la salida del plan se registran en esta validación; los artifacts locales completos se mantuvieron temporales para no añadir una segunda carpeta permanente de evidencias.

## Pendiente / estado de lanzamiento

- La suite monolítica no se repitió: su baseline agotó el límite externo de 480 s sin resumen, así que no se declara ese comando PASS. La matriz oficial de shards sí está verde localmente.
- Actions anterior `39689871971` corresponde al SHA previo `aff565a`: shard 6 abortó; `windows-smoke` y `flatpak-smoke` pasaron. Debe verificarse el run del SHA que contenga este fix antes de cerrar el gate remoto.
- La auditoría de los cuatro paquetes y la evidencia de aceptación manual siguen pendientes. `docs/release/manual-acceptance-0.3.0.md` se conserva y los casos sin evidencia continúan `NOT TESTED`.
- No merge, tag, Release ni publicación. Estado de preparación beta permanece **NO-GO** hasta completar los gates de distribución/aceptación y revisar el Actions run del SHA final.
