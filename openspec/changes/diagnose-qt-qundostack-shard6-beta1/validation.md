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

## Actions sobre el SHA final del fix

Run [#39697199469](https://github.com/PJGV333/Chemuson/actions/runs/39697199469), SHA `8c28d47c0518f4ca40619839ae3504df8cf87096`:

- `pytest-plan` PASS (2.598 IDs); `windows-smoke` PASS (33 s); `flatpak-smoke` PASS (25 s).
- Shards 0–4 y 7 PASS. Shard 5 agotó el límite externo de 300 s: 158/268 timings, todos pasados hasta `test_build_appimage_preview_creates_type2_without_public_updater_metadata`; JUnit ausente.
- Shard 6 ya no abortó en onboarding: `test_onboarding_close_releases_overlay` y el resto de `test_ui_polish.py` pasaron. Sin embargo, también agotó 300 s: 342/373 timings (340 passed, 2 skipped), último node iniciado `tests/test_ui_svg_icons.py::TestDynamicBuilders::test_atom_escapes_special_chars`; JUnit ausente.
- `pytest-summary` rechazó la colección por esos dos shards incompletos. No se aumentó el timeout ni se movieron tests. No se repitieron esos comandos exactos tras timeout; cobertura completa ya está comprobada en la matriz local secuencial de 2.598 IDs. Los artifacts remotos de los shards tienen el log/timings/process report.

## Estado de lanzamiento

El fix del SIGABRT está validado por red/green A/B, por la matriz local completa y por el hecho de que el test de onboarding pasa en Actions. **La verificación oficial remota aún NO está verde** por dos timeouts bajo el límite existente, aunque Windows y Flatpak smoke sí pasan. No declarar CI completo PASS.

La suite monolítica no se repitió: su baseline agotó el límite externo de 480 s sin resumen. La auditoría de los cuatro paquetes y la evidencia de aceptación manual siguen pendientes. `docs/release/manual-acceptance-0.3.0.md` se conserva y los casos sin evidencia continúan `NOT TESTED`.

Estado de preparación beta: **NO-GO** hasta resolver mediante alcance/decisión autorizada los timeouts de Actions y completar distribución/aceptación. No merge, tag, Release ni publicación.

## Addendum — SIGABRT de shard 5 y perfil actualizado (2026-10-21)

### Red/green ordenado

Antes de editar, el estado Git fue `stabilize/v0.3.0-beta.1` en `3f800c715ee937d9ac72385d8189bb9d96902ea3`; los tres cambios pendientes del perfil/contratos se conservaron y su diff se respaldó en `/tmp/chemuson-beta1-pre-diagnostic-3f800c7.patch` (SHA-256 `47a8a5e1f3aadcfe3735034ecf0aec9501d183cba74a0c7899f10a56d6019083`). El área `.tmp-beta1-readiness/` se inspeccionó y retuvo durante el diagnóstico.

Se ejecutó exactamente la asignación/orden de 349 IDs del shard 5 del plan candidato en snapshot base de `git archive` de `3f800c7`. Base sin fix: **SIGABRT/exit -6**, 272 timings, JUnit ausente, 93.597 s; último node iniciado `tests/test_molecular_assistant_lifecycle.py::test_unstarted_dialog_close_and_window_close_leave_no_registry_entries`. El stack nativo atraviesa `_process_qt_deletes`, `QUndoStack::~QUndoStack` y `cleanChanged`. `chemuson.__file__` quedó confirmado dentro del snapshot.

La variante corrigió únicamente el receiver de limpieza en tests. El mismo snapshot/orden completó **341 passed, 8 skipped**, 349/349 JUnit y timings, exit 0 en 151.017 s. El archivo lifecycle completó 9/9 en 3.95 s. En el checkout el shard 5 repitió 349/349 en 151.145 s; tras limitar también los probes equivalentes, los focales fueron: lifecycle + `test_ui_polish.py` **47 passed**, interoperabilidad MOL + shutdown subprocess **16 passed**.

Causa confirmada: `_process_qt_deletes()` vaciaba globalmente `DeferredDelete` de objetos Qt ajenos en un proceso que comparte QApplication. La función final requiere el receiver, envía el evento directamente antes de `processEvents()` y así conserva las aserciones existentes de destrucción, registros vacíos, cierre y workers. El primer experimento que llamaba `processEvents()` antes del envío scoped dio `RuntimeError` porque Qt ya había destruido el diálogo; se corrigió el orden del helper, no se debilitó ninguna aserción. Las llamadas de `tests/probes/p0_shutdown_workflow.py` y `tests/probes/p0_open_mol.py` se inspeccionaron y también apuntan al objeto específico; sus pruebas corrieron en subprocess dedicados. No se modificó producción.

### Validación final: perfil candidato + Qt

Perfil validado: `tests/ci/pytest_file_durations.json`, SHA-256 `d2d7db19a9e651d59170d2548dadd8b551623ab4468d99391a29939feea5341e`; procede de Actions #39697199469 (`8c28d47c`), con #39679607676 y #38996089270 como suplementarios. Tiene cobertura exacta de los 253 archivos del manifiesto, fallback 0.2 s/test. Los contratos del perfil/workflow: **10 passed**. Plan reproducible: 2.598 IDs, distribución `[257,371,353,261,316,349,296,395]`, carga estimada ≈207.84 s/shard.

Se ejecutaron secuencialmente los ocho shards con límite del runner existente de 300 s: `[72.072,144.049,84.579,81.085,100.477,152.034,65.211,177.566]` s. Resultado **2.576 passed, 22 skipped, 0 failed, 0 errors**; JUnit y timings 2.598/2.598; verifier confirma ocho shards completos. El máximo fue 177.566 s. El manifiesto/hash `b2cb90fd…7248f` no cambió. No hubo skips/xfails nuevos ni expansión de timeout.

Los archivos de salida completos se preservaron fuera del repositorio bajo `/tmp/chemuson-beta1-readiness-3f800c7/`. El manifiesto raíz `EVIDENCE-MANIFEST.sha256` cubre 360 archivos (SHA-256 `61046269db7c8bd4f53a648eebd809f84deebdc96d62a2ac5a93d0a333f3475a`). Los snapshots temporales fuente se retiraron después de preservar los logs/reportes.

## Addendum — run Actions posterior al receiver-scoping (2026-10-22)

El run histórico descrito arriba (`39697199469`) queda supersedido para el SHA final por [Actions #39723284720](https://github.com/PJGV333/Chemuson/actions/runs/39723284720), push de `stabilize/v0.3.0-beta.1` sobre `5ca26395a125acbabb38573d38c0605504b57cb2`. Resultado global: **success**; `pytest-plan`, ocho shards, `pytest-summary`, `windows-smoke` y `flatpak-smoke` pasaron.

El plan/JUnit/timings cubrieron exactamente **2.598/2.598 IDs**: **2.576 passed, 22 skipped, 0 failed/errors/timeouts**. El shard 5, cuya secuencia local había expuesto el SIGABRT en el cierre de Molecular Assistant, completó sus 349 asignaciones en 247.689 s: 341 passed y 8 skipped, JUnit/timings completos, exit 0. El shard 7 fue el más lento con 293.063 s, dentro del timeout existente de 300 s. No se aumentó el límite ni se alteraron aserciones o skips.

Esto cierra la validación remota del diagnóstico/fix de los helpers Qt y del perfil CI para este SHA. El margen del shard 7 es de sólo 6.937 s y se conserva como riesgo operativo. **No es aceptación manual ni aprobación de publicación**: siguen pendientes auditoría de los cuatro formatos y ejecución de `docs/release/manual-acceptance-0.3.0.md`. El dictamen de beta continúa **NO-GO**; sin merge, tag, Release o publicación.

## Addendum — follow-up CI (2026-10-22)

Aunque Actions #39723284720 sobre SHA `5ca26395` completó la matriz y ambos smoke jobs, el follow-up #39724355063 sobre `a67c9152` terminó incompleto por timeouts en shards 5 y 7 (302.047 s / 302.071 s, SIGTERM, JUnit ausente en ambos). Entre los IDs que sí terminaron no se reportaron fallos; el summary job rechazó correctamente la cobertura incompleta. La ejecución previa sigue siendo evidencia verde para el SHA de código `5ca26395`, pero no se debe extrapolar como PASS del HEAD remoto actual. Estos timeouts no reprodujeron el SIGABRT `QUndoStack` que corrige este cambio y no se les atribuye esa causa. El receiver-scoping permanece validado; el diagnóstico de shard 5 continúa cerrado, mientras CI general/perfil quedan bloqueados para investigación separada. **NO-GO**; no se aumentó timeout, no se hizo merge, tag ni Release.
