# Validación — estabilización final ChemUSON 0.3.0-beta.1

## Cambios y diagnóstico

- Rama local `stabilize/v0.3.0-beta.1`, base `61f7647b683db3eca19ea06e5ec82ab2897d640a`. Sin commit, push, tag, Release, PR, worktree ni publicación. Cambios de código limitados al controller CompChem y sus regresiones; nuevo node ID registrado explícitamente en `tests/ci/expected_pytest_nodeids.txt`.
- `HeldCompChemWorker` ahora declara `progress_ready`, señal que `CompChem3DController.start_job` conecta. El test de cierre deja de fallar antes de sus aserciones.
- El fallo CI CompChem tenía `worker.finished` observado, `pending_results=(1,)` y el QThread aún activo al límite preexistente de 3 s. La conexión normal `worker.finished -> thread.quit` encolaba `quit` hacia el hilo principal. Regresión roja previa al fix: el worker emite `finished`, pero un main thread bloqueado en `QThread.wait(1000)` no puede entregar el quit; tras bombear eventos, el thread termina. Fix mínimo: conectar sólo `thread.quit` con `Qt.ConnectionType.DirectConnection` (thread-safe); callbacks de resultado/cierre siguen en el event loop GUI. No se amplió el timeout.
- La regresión nueva pasa en Python 3.11.15/PyQt6 6.11.0/Qt runtime 6.11.2 y en el Python del repo. Esto identifica y corrige un modo de fallo compatible con el estado CI; **no sustituye un nuevo Actions run en este SHA**.

## Verificación focal y estática

- Python del repo: Python 3.14.8, PyQt6 6.11.0/Qt 6.11.2, RDKit 2026.03.6, pytest 9.1.1.
- Venv temporal `/tmp/chemuson-final-stabilization-311`: Python 3.11.15, PyQt6 6.11.0/Qt 6.11.2, RDKit 2026.9.1, pytest 9.1.1; coincide en versiones de paquetes con el runner CI #39679607676 salvo patch de Python.
- Rojo/verde del nuevo test `test_compchem_worker_stops_its_thread_without_main_event_loop`: **1 failed** antes del cambio de producción (`thread.wait(1000)==False`, `worker.finished` sí observado); **1 passed** tras el `DirectConnection`.
- Python 3.11, tests CompChem + regresión: **2 passed**. Tests focales del doble de cierre: **1 passed**. Wrapper subprocess `test_existing_async_worker_shutdown_contract_in_isolated_pytest_process`: **1 passed**. Ruff focal: PASS.
- Python del repo, `tests/test_compchem3d_dock.py tests/test_gui_async_worker_shutdown.py tests/test_p0_qt_shutdown_subprocess.py`: **62 passed, 1 skipped** en 63.68 s. Python del repo `tests/test_ui_polish.py` aislado: **38 passed** en 11.20 s. `tests/test_test_workflow_contract.py`: **4 passed**. `tests/architecture`: **283 passed**.
- `python -m compileall -q src tests tools packaging`: PASS. `openspec validate final-stabilization-chemuson-0-3-0-beta-1 --strict`: PASS. `git diff --check`: PASS. Ruff global continúa con sólo el F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se tocó.
- Se actualizó el manifiesto de node IDs únicamente tras revisar diff exacto: 2.597 → **2.598**, adición única `tests/test_compchem3d_dock.py::test_compchem_worker_stops_its_thread_without_main_event_loop`, sin eliminaciones. `test_test_workflow_contract.py` y el plan confirman el manifiesto.

## Matriz local de ocho shards

- Plan existente `tools/ci_pytest_shards.py`, ejecución serial, Python del repo, timeout externo 300 s/shard. Plan: 2.598 IDs; distribución `[248,243,327,362,449,268,373,328]`; hash del manifiesto `b2cb90fdc74a1c2417885557f9391f12924235c3043d2e1a73fee0340027248f`. Evidencia en `/tmp/chemuson-final-stabilization-shards/`.
- Shards 0, 1, 2, 3, 4, 5 y 7 terminaron exit 0 con JUnit/timings completos: conjuntamente **2.205 passed, 20 skipped**.
- Shard 6: **SIGABRT/exit -6**, no timeout, a los 144.659 s. JUnit no se escribió; timings registran 308/373 IDs (**306 passed, 2 skipped**) y 65 IDs sin registro completo. Último test iniciado: `tests/test_ui_polish.py::test_onboarding_close_releases_overlay`. La pila nativa acaba en `QUndoStack::~QUndoStack` / `QUndoStack::cleanChanged` durante destrucción de QWidget. El shard se repitió una vez en directorio separado: mismo SIGABRT en el mismo test (143.567 s, mismos 308 timings). El archivo de UI aislado sí pasa 38/38; por tanto el problema depende del orden/estado Qt previo. No se atribuye al cambio CompChem sin prueba ni se amplía el alcance con una refactorización de UI.
- `tools/ci_pytest_shards.py verify` rechazó correctamente la matriz: job status `failure`, shard 6 exit -6/SIGABRT y JUnit ausente. **No es una matriz verde**. Los siete shards verdes y el shard abortado cubren el intento, no los 2.598 IDs completos.
- Un test local adicional, `test_file_workflow_close_and_autosave_recovery_in_bounded_processes`, falló en la venv Python 3.11 porque su helper añade primero `chem/lib/python3.14/site-packages` (ruta hardcodeada del repo) al `PYTHONPATH`; RDKit falla con ABI `undefined symbol: PyObject_DelAttr`, seguido de abort del subproceso. El mismo test pasa con el Python 3.14 del repo (**1 passed**). No se modificó el helper ni se trata como defecto de producción; se registra como incompatibilidad del runner Python3.11 temporal local.

## Shards del Actions baseline y skips

El run remoto existente #39679607676 sigue siendo del SHA base `61f7647b683db3eca19ea06e5ec82ab2897d640a`: **2.572 passed, 22 skipped, 3 failed**. Los dos fallos de doble/CompChem se corrigen localmente como arriba, pero no se ha disparado un run nuevo y el resultado remoto sigue `failure`.

Los 22 skips JUnit del Actions baseline, con su razón literal:

1. `tests/test_block_unwrap_depiction.py::test_tetrandrine_import_generates_block_unwrap_candidate` — RDKit depiction was not donut-like enough to trigger block unwrap.
2. `tests/test_compchem3d_dock.py::test_promoted_explicit_hydrogen_is_not_duplicated_by_rdkit` — RDKit no disponible en este entorno.
3. `tests/test_scaffold_depiction.py::test_tetrandrine_best_depiction_not_forced_to_single_strategy` — RDKit tetrandrine depiction did not trigger scaffold candidate.
4. `tests/test_scaffold_depiction.py::test_vancomycin_best_depiction_has_scaffold_candidate` — RDKit vancomycin depiction did not trigger scaffold candidate.
5–7. `tests/test_rdkit_strict_validation.py::{test_coordination_bond_is_exported_as_dative,test_strict_validation_accepts_tetramethylammonium,test_strict_validation_rejects_overvalent_neutral_carbon}` — RDKit no disponible.
8–9. `tests/test_chemname_pr27.py::{test_stereochemistry_e_2_butene,test_stereochemistry_lactic_acid}` — RDKit no disponible.
10. `tests/test_geometry3d_service.py::test_obabel_backend_degrades_or_runs` — Open Babel no disponible.
11–15. `tests/test_chemname_acceptance.py::test_acceptance_case[smiles_triazole_123]`, `[smiles_triazole_124]`, `[smiles_tetrazole]`, `[smiles_benzotriazole]`, `[smiles_nickel_tetracarbonyl]` — `worker_error:invalid_input`.
16–19. `tests/test_rdkit_roundtrip.py::{test_duplicate_bonds_do_not_crash_rdkit_conversion,test_molgraph_roundtrip_smiles,test_pseudoatoms_use_dummy_atoms_in_rdkit_conversion,test_smiles_import_regression_for_naphthyl_urea}` — RDKit no disponible.
20–21. `tests/test_stereo_advanced.py::{test_allene_axial_best_effort,test_biaryl_atropisomeric_best_effort}` — RDKit no disponible.
22. `tests/test_viewport3d_depth_occlusion.py::test_real_viewport_occludes_bonds_with_depth_after_qpainter_overlays` — requires an authorized Qt/OpenGL desktop session.

La ejecución local Python 3.14 tuvo 20 skips: coincide en estos motivos salvo que los dos casos de depiction de tetrandrina y explicit hydrogen se ejecutaron en vez de skip con esa versión de dependencias. Los skips no son aceptación manual.

## Distribución, aceptación manual y dictamen

- **No auditado aún en esta continuación:** cuatro builds nuevos (Windows portable, Inno Setup, AppImage Type 2, Flatpak), provenance/checksums de un mismo SHA, ChemName/RDKit/3D congelados, arranque/extracción/FUSE y aislamiento Flatpak. No existe un preview nuevo para el SHA estabilizado; todos permanecen **NOT BUILT / NOT TESTED** aquí.
- `docs/release/manual-acceptance-0.3.0.md` no se alteró ni se convirtió en PASS. Sin artifacts identificados y observación manual autorizada, los casos gráficos/químicos y Windows/OpenGL permanecen según su estado previo (NOT TESTED donde no hay evidencia).
- **Dictamen: NO-GO / bloqueado.** CI remoto continúa rojo en el SHA base; matriz local actual aborta reproduciblemente en Qt teardown; faltan run remoto y cuatro builds exactos, además de aceptación manual. No cambiar estos estados sin decisión del propietario sobre el aborto `QUndoStack` y evidencia nueva.
- Ningún push, commit, tag, Release, PR, preview, actualización de canal ni cambio a `main`/`gh-pages` fue realizado.

## Addendum — Qt shard 5 + perfil medido (2026-10-21)

Se diagnosticó el SIGABRT que introdujo el perfil candidato al reasignar lifecycle a shard 5. Red/green del mismo plan de 349 IDs sobre snapshot base `3f800c7`: antes, exit -6 tras 272 timings y JUnit ausente en `test_unstarted_dialog_close_and_window_close_leave_no_registry_entries`; después de limitar `DeferredDelete` al QObject exacto, 341 passed + 8 skipped, JUnit/timings completos y exit 0. La causa fue el flush global `sendPostedEvents(None, DeferredDelete)` en el helper de lifecycle; no se tocó producción. También se inspeccionaron y acotaron los probes Qt equivalentes, preservando todas las aserciones.

El perfil local se mantiene: fue reconstruido desde métricas de Actions documentadas, cubre exactamente 253 archivos y usa fallback 0.2 s/test. Los contratos pasan **10/10**. El plan de ocho shards balancea carga estimada en ~207.84 s por shard. Matriz secuencial final: **2.598/2.598 IDs, 2.576 passed, 22 skipped, 0 failed/errors**, JUnit/timings exactos; wall time individual máximo **177.566 s < 300 s**. Se preservan orden, manifiesto, timeout y aserciones.

Este estado se registró antes del push y queda supersedido respecto al CI por el addendum del 2026-10-22. Los builds de los cuatro paquetes y la aceptación manual continúan pendientes; la publicación permanece **NO-GO**.

## Addendum — Actions remoto de la revisión final (2026-10-22)

Se publicaron mediante push normal, sólo a `stabilize/v0.3.0-beta.1`, los commits `9ced6c341f673e7e1512d6ba9e1afb97f327543c` (receiver-scoping de `DeferredDelete` en tests/probes Qt) y `5ca26395a125acbabb38573d38c0605504b57cb2` (perfil CI derivado de timings y contratos). El SHA remoto de la rama quedó en `5ca26395`; `main` y `release/v0.3.0-beta.1-prep` no se modificaron.

Actions [#39723284720](https://github.com/PJGV333/Chemuson/actions/runs/39723284720), ejecutado sobre el SHA exacto `5ca26395a125acbabb38573d38c0605504b57cb2`, terminó **success**. `pytest-plan`, los ocho shards, `pytest-summary`, `windows-smoke` y `flatpak-smoke` aprobaron. El plan conservó 2.598 IDs y el hash esperado `b2cb90fdc74a1c2417885557f9391f12924235c3043d2e1a73fee0340027248f`; cada asignación coincide con JUnit y timings y se ejecutó exactamente una vez.

| Shard | IDs/JUnit | Passed | Skipped | Tiempo de pytest |
|---:|---:|---:|---:|---:|
| 0 | 257/257 | 256 | 1 | 114.959 s |
| 1 | 371/371 | 371 | 0 | 198.091 s |
| 2 | 353/353 | 349 | 4 | 186.283 s |
| 3 | 261/261 | 261 | 0 | 99.836 s |
| 4 | 316/316 | 310 | 6 | 167.885 s |
| 5 | 349/349 | 341 | 8 | 247.689 s |
| 6 | 296/296 | 293 | 3 | 117.191 s |
| 7 | 395/395 | 395 | 0 | 293.063 s |

Total: **2.576 passed, 22 skipped, 0 failed, 0 errors**, JUnit y timings **2.598/2.598**, ningún timeout ni señal. El shard más lento quedó 6.937 s por debajo del límite existente de 300 s. Aunque CI aprobó, el margen es estrecho; se conserva el límite y el perfil sin cambios y se registra como riesgo de capacidad para la siguiente ejecución, no como motivo para ocultar/relajar pruebas.

Los artifacts remotos descargados y la evidencia local preservada (360 archivos de evidencia, más el manifiesto actualizado tras el follow-up) están fuera del worktree en `/tmp/chemuson-beta1-readiness-3f800c7/`; `EVIDENCE-MANIFEST.sha256` verifica 360/360 archivos (SHA-256 del manifiesto `61046269db7c8bd4f53a648eebd809f84deebdc96d62a2ac5a93d0a333f3475a`). Los snapshots temporales fuente se retiraron tras preservar los logs y reportes. El diff previo se conserva en `/tmp/chemuson-beta1-pre-diagnostic-3f800c7.patch` (SHA-256 `47a8a5e1f3aadcfe3735034ecf0aec9501d183cba74a0c7899f10a56d6019083`).

Este PASS acredita CI automatizado del SHA de código publicado, no auditoría/instalación de cuatro paquetes ni aceptación manual. `docs/release/manual-acceptance-0.3.0.md` permanece intacto; los casos sin observación siguen `NOT TESTED`. Por tanto, el dictamen de beta permanece **NO-GO para publicación** hasta completar la auditoría de distribución y la aceptación pendiente; no hubo merge, PR, tag, Release ni publicación.

## Addendum — follow-up de Actions y prueba de rebalanceo (2026-10-22)

El commit documental `a67c915275d407ff0a8081f1f7ee1a8ef9dbcaa2` se publicó normalmente después de registrar el run #39723284720; no cambió código, tests ni perfil respecto a `5ca26395`. Sin embargo, su Actions [#39724355063](https://github.com/PJGV333/Chemuson/actions/runs/39724355063) terminó **failure**: `pytest-plan`, Windows smoke y Flatpak smoke pasaron, pero shards 5 y 7 recibieron SIGTERM por el límite existente de 300 s; el resumen fail-closed rechazó la matriz. No hubo fallos pytest registrados entre las pruebas completadas, pero el resultado es incompleto y no es PASS.

| Shard | Asignadas | Timings completos | Tiempo | Último test iniciado | Resultado |
|---:|---:|---:|---:|---|---|
| 5 | 349 | 346 | 302.047 s | `tests/test_validation_dock.py::test_validation_controller_reduces_selected_bond_when_safe` | timeout/SIGTERM, sin JUnit; faltan 3 IDs |
| 7 | 395 | 335 | 302.071 s | `tests/test_sp3_bond_placement.py::test_default_bond_endpoint_ignores_isolated_hidden_carbon_placeholder` | timeout/SIGTERM, sin JUnit; faltan 60 IDs |

El run conservó 2.535 timing records (2.513 passed, 22 skipped; 0 failures observadas), no 2.598 resultados completos. Es la misma asignación/perfil que pasó en #39723284720, lo que evidencia sensibilidad de duración; no se atribuye a un test colgado ni se amplió el límite. El SHA de código `5ca26395` sí tiene una ejecución completa verde (#39723284720), pero el HEAD remoto actual `a67c9152` sólo tiene este follow-up incompleto.

### Experimento de perfil no adoptado

A partir de los promedios por archivo del run completo #39723284720 y los timings finalizados del run parcial #39724355063 se generó, sin instalarlo como cambio aceptado, un candidato de perfil SHA-256 `fdb58a09032176c33ada647cc4200ff092bb41b80e8b205807b3f7cfdfb3254e`. Conserva fallback 0.2 s/test y cobertura de 253 archivos; reasigna 2.598 IDs en `[249,280,388,309,373,357,386,256]`, con cargas estimadas de 268.72–268.74 s. Contratos de perfil/workflow: **10 passed**; ciclo de vida/UI focal **47 passed**; wrappers P0 **16 passed**.

La matriz local con ese perfil experimental también fue incompleta: shards 0,1,3–7 pasaron; shard 2 agotó 300 s (302.018 s, 365/388 timings, exit -15/SIGTERM, sin JUnit), con último node iniciado `tests/test_ui_theme_foundation.py::TestApplyTheme::test_apply_theme_on_application`. Ese nodo pasó aislado (1/1) y el archivo completo pasó (28/28), así que la causa de la condición orden-dependiente no queda clasificada. No se repitió el shard exacto ni se movió el test para ocultar el resultado. El candidato no se conservó en el árbol de trabajo; queda archivado junto con plan, logs, timing files, lista de IDs incompletos y generador en `/tmp/chemuson-beta1-readiness-3f800c7/profile-rebalance-baseline/` y `local-rebalanced-profile/`.

**Estado:** no aceptar el perfil experimental ni declarar CI verde para HEAD `a67c9152`. Se requiere una decisión para diagnosticar por separado el orden Qt que acaba en `test_apply_theme_on_application` antes de otra redistribución; no tocar producción, aserciones, skips/xfails ni el límite de 300 s. Publicación beta continúa **NO-GO**; no se hizo merge, PR, tag, Release ni publicación.

## Addendum — registro oficial CI #422 y auditoría de distribución (2026-10-22)

### CI del SHA estabilizado

El run que el owner identifica como **CI #422** corresponde al GitHub Actions run [#39782991976](https://github.com/PJGV333/Chemuson/actions/runs/39782991976), workflow `test`, evento `push`, rama `stabilize/v0.3.0-beta.1`, SHA exacto `f0151cf09a7200acc53cb6b371105e79fd3cac1d`. Conclusión `success`; duración de pared 6:25. `pytest-plan`, ocho shards, `pytest-summary`, `windows-smoke` y `flatpak-smoke` terminaron correctamente. Estos dos smoke jobs no compilan ni instalan el set de cuatro paquetes.

El plan asignó 2.598 node IDs una sola vez, hash del manifiesto `b2cb90fdc74a1c2417885557f9391f12924235c3043d2e1a73fee0340027248f`, perfil `d2d7db19a9e651d59170d2548dadd8b551623ab4468d99391a29939feea5341e`. Los JSON por shard, JUnit y timings se descargaron; cada JUnit/timing coincide con su asignación, sin duplicados ni IDs ausentes.

| Shard | IDs/JUnit/timings | Passed | Skipped | Failed/errors | Tiempo pytest |
|---:|---:|---:|---:|---:|---:|
| 0 | 257/257/257 | 256 | 1 | 0/0 | 99.329 s |
| 1 | 371/371/371 | 371 | 0 | 0/0 | 299.398 s |
| 2 | 353/353/353 | 349 | 4 | 0/0 | 140.088 s |
| 3 | 261/261/261 | 261 | 0 | 0/0 | 161.943 s |
| 4 | 316/316/316 | 310 | 6 | 0/0 | 189.076 s |
| 5 | 349/349/349 | 341 | 8 | 0/0 | 209.207 s |
| 6 | 296/296/296 | 293 | 3 | 0/0 | 75.786 s |
| 7 | 395/395/395 | 395 | 0 | 0/0 | 289.271 s |
| **Total** | **2.598/2.598/2.598** | **2.576** | **22** | **0/0** | — |

El shard 1 quedó a **0.602 s** del límite externo de 300 s. Esto es una ejecución completa válida, pero no una demostración estadística de estabilidad. El run anterior [#39724355063](https://github.com/PJGV333/Chemuson/actions/runs/39724355063) sobre `a67c9152` —SHA previo cuya diferencia posterior a este run es documental— tuvo timeout en shards 5 y 7 a ~302 s. No se aumentó el límite, no se aceptó el perfil experimental `fdb58a09…` y no se omitieron pruebas.

La matriz automatizada incluyó, entre otros: `tests/test_amide_urea_2d_continuity.py` 19/19; MOL P0 14/14; importación estéreo SMILES 12/12; dibujo en bordes 6/6; `drag_move_undo` 3/3; Live MM stream 13/13; proyección estereoquímica 5/5; CompChem 3D 58/59 (un skip RDKit); servicio geometry3D 35/36 (un skip Open Babel); contratos de worker RDKit empaquetado 8/8 (no ejecutan los artifacts del preview); SVG 41/41; tema 28/28; ChemName acceptance 71/76. `clean2d_tetrandrine_like_local` pasó 5/5; los smokes candidatos de block-unwrap y scaffold tuvieron 1 y 2 skips geométricos respectivamente. El run también ejercitó el corpus general Clean2D, edición 2D, guardado, undo/redo, 3D y tests de UI dentro de los otros shards. Estos números describen tests fuente/automatizados en pytest, no revisión química independiente ni aceptación gráfica manual.

Los 22 skips JUnit del run (razón literal resumida) fueron:

- Qt/OpenGL no disponible: `tests/test_viewport3d_depth_occlusion.py::test_real_viewport_occludes_bonds_with_depth_after_qpainter_overlays` (requiere sesión autorizada Qt/OpenGL).
- Geometría heurística no activada: `tests/test_block_unwrap_depiction.py::test_tetrandrine_import_generates_block_unwrap_candidate` (no se detectó donut); `tests/test_scaffold_depiction.py::test_tetrandrine_best_depiction_not_forced_to_single_strategy` y `::test_vancomycin_best_depiction_has_scaffold_candidate` (RDKit no produjo candidato scaffold).
- RDKit no disponible en los tests de entorno: `tests/test_rdkit_strict_validation.py::{test_coordination_bond_is_exported_as_dative,test_strict_validation_accepts_tetramethylammonium,test_strict_validation_rejects_overvalent_neutral_carbon}`; `tests/test_chemname_pr27.py::{test_stereochemistry_e_2_butene,test_stereochemistry_lactic_acid}`; `tests/test_rdkit_roundtrip.py::{test_duplicate_bonds_do_not_crash_rdkit_conversion,test_molgraph_roundtrip_smiles,test_pseudoatoms_use_dummy_atoms_in_rdkit_conversion,test_smiles_import_regression_for_naphthyl_urea}`; `tests/test_compchem3d_dock.py::test_promoted_explicit_hydrogen_is_not_duplicated_by_rdkit`; `tests/test_stereo_advanced.py::{test_allene_axial_best_effort,test_biaryl_atropisomeric_best_effort}`.
- Casos ChemName con `worker_error:invalid_input`: `tests/test_chemname_acceptance.py::test_acceptance_case[smiles_triazole_123]`, `[smiles_triazole_124]`, `[smiles_tetrazole]`, `[smiles_benzotriazole]`, `[smiles_nickel_tetracarbonyl]`.
- Open Babel opcional ausente: `tests/test_geometry3d_service.py::test_obabel_backend_degrades_or_runs`.

Ninguno de esos skips se cambió por esta tarea. La aceptación P0/P1 manual, OpenGL real y ejecución en Windows instalado siguen separadas.

### Cuatro artefactos Build Preview existentes — SHA no candidato

Se descargó e inspeccionó el Build Preview [#39679607716](https://github.com/PJGV333/Chemuson/actions/runs/39679607716). Su reporte es `complete=true`, los cuatro jobs `success`, `publication=false`, versión `0.3.0-beta.1`, pero **el source SHA es `61f7647b683db3eca19ea06e5ec82ab2897d640a`**, de `release/v0.3.0-beta.1-prep`, no el SHA de estabilización `f0151cf…`. Artifact IDs, tamaños y hashes SHA-256 calculados localmente:

| Formato | Artifact Actions | Archivo | Bytes | SHA-256 del paquete | SHA-256 de provenance | SHA-256 del ZIP Actions |
|---|---:|---|---:|---|---|---|
| Windows portable | [12383692156](https://github.com/PJGV333/Chemuson/actions/runs/39679607716/artifacts/12383692156) | `Chemuson-v0.3.0-beta.1-preview-61f7647b-windows-x86_64-portable.exe` | 133,250,798 | `178e27bf444ef088a9386a79c36349640c7fb1c7f96b4d7451d4024f5a1169a2` | `a18dfeee9c30ab398e41e76db644fdd910a59008fe650010bd9364b691310d9f` | `b38501e8dc25a3d15e111c164938f603c3e2c97f0cf9cbaa2d994454a8c07dd2` |
| Windows Inno Setup | [12383775788](https://github.com/PJGV333/Chemuson/actions/runs/39679607716/artifacts/12383775788) | `Chemuson-v0.3.0-beta.1-preview-61f7647b-windows-x86_64-setup.exe` | 134,386,949 | `cfedeaf05d65c9a27193736f80562ec2a50c15e86bcfaef1404a89f5a344c113` | `aedbab3d5dabdf04096fa0bb0767126b758e28d82cbe95af56d896400d87e624` | `460a5b514bb12451d9dcdd6e2a3f8280c7a73a9a141bcbf112ad6679e0036af5` |
| Linux AppImage Type 2 | [12384370204](https://github.com/PJGV333/Chemuson/actions/runs/39679607716/artifacts/12384370204) | `Chemuson-v0.3.0-beta.1-preview-61f7647b-linux-x86_64.AppImage` | 176,645,312 | `07f5cb0639336c0c120fc2671537ff645db66bbb49d735a030731182379e7a0f` | `6ce5a4b3075a4450a76e3437670e34f842a843b8bf9dbc4ce9b12be8d3207c0c` | `9ccbffe659314ee1d419662eef0a1033c96048f383e99d2c5b97f062e9ca7f3b` |
| Linux Flatpak | [12384500783](https://github.com/PJGV333/Chemuson/actions/runs/39679607716/artifacts/12384500783) | `Chemuson-v0.3.0-beta.1-preview-61f7647b-linux-x86_64.flatpak` | 110,085,312 | `4fe0d0a94989f45cb92751651d3b12d0ccd9251a5ce17ad4d4876fdf526b7f26` | `1d2a1f2c4f2b778058ee5c5721d50c432305f33e613b6eb5b6fc11020422729c` | `c602081942d9eb3a3825fcda8e1bde03baed7c63d17afdbe8e739905eab89f3f` |

Los cuatro hashes internos y tamaños coinciden con `checksums.sha256` y `preview-provenance.json` de su grupo. El reporte está en artifact [12384141020](https://github.com/PJGV333/Chemuson/actions/runs/39679607716/artifacts/12384141020). Las ZIP SHA-256 se tomaron del log de subida Actions; los paquetes no son un set para el SHA estabilizado actual.

**Windows portable — inspección estática y smoke automatizado:** el PE es x86-64 GUI. El runner Windows ejecutó icon smoke en tema claro/oscuro y DPR 2: 69/69 SVG con píxeles visibles y 11 controles esenciales. El smoke congelado RDKit verificó `rdchem`, `rdMolDescriptors` y `rdDistGeom` dentro de `_MEIPASS`, logP -0.0014, TPSA 20.23, HBD/HBA 1, SMILES `CCO` y conformero 3D; parent sin módulos RDKit/GUI. ChemName comprobó 9/9 templates y 6/6 nombres frente a Python. No es una sesión gráfica/manual del propietario.

**Windows Inno Setup — build, no instalación:** el log registra compilación Inno Setup 6.7.1 exitosa en 21.078 s, incluyendo el `dist/Chemuson.exe` validado y marker; setup tiene stub PE GUI x86 (payload portable x64). No se instaló, actualizó, desinstaló ni se inspeccionó en Windows interactivo. Las pruebas locales `7z` no reconocen el contenedor Inno; esto no invalida el build del runner, pero deja inspección/ejecución nativa manual pendiente.

**Linux AppImage — Type 2 y bundle smoke automatizados:** inspección local confirmó ELF64 x86-64 y marker `AI\x02`; el paquete se extrajo mediante `--appimage-extract` sin FUSE, `AppRun --version` devolvió `0.3.0-beta.1`, y los sidecars públicos `.updateinfo`, `.update.json`, `.zsync` están ausentes como corresponde al preview. El job remoto ejecutó `validate_appimage.py` completo con éxito: AppDir, desktop/AppStream, PyInstaller resources, 69 SVG, ChemName 9/9 + 6/6, RDKit/extensiones/descriptores/SMILES/3D y headless launch. El intento local del validator terminó en `Could not inspect the bundled PyInstaller CArchive` porque el intérprete anfitrión carece del lector PyInstaller; no se altera el resultado remoto y no se declara validación CArchive local.

**Linux Flatpak — import/recursos, sin instalación en la sesión:** el log de Build Preview construyó el bundle y ejecutó ChemName smoke con `/app/lib/python3.13/site-packages`: 9/9 recursos y 6/6 controles. Se importó además el bundle descargado a un repositorio OSTree temporal local; ref `app/io.github.PJGV333.Chemuson/x86_64/preview-61f7647b`, commit `afb42bae27f70875651ad11abb09055e08c32ac3b8a3ca757e2a02853a3a002d`, y checkout en modo user encontró los nueve `.mol`. No se instaló en perfil usuario/sistema, no se añadió remoto ni se inició GUI Flatpak.

### Identidad del candidato, límites y dictamen

El Actions run #39782991976 es el test workflow y **no generó paquetes**. El último Build Preview descargado (#39679607716) corresponde a SHA `61f7647b…`. `git diff 61f7647b… f0151cf… -- src/chemuson packaging chemuson.spec` muestra un cambio de producción en `src/chemuson/gui/controllers/compchem3d_controller.py` (fix de salida de QThread); por tanto, esos binarios no pueden presentarse como reproducibles desde el nuevo SHA.

`build-preview.yml` sólo acepta push de `release/**-prep` o dispatch cuya validación requiere ese patrón. El HEAD actual está en `stabilize/v0.3.0-beta.1`. No se creó una rama puente, no se cambió el filtro, no se usó un ref ajeno y no se lanzó un preview inválido. Para construir cuatro artefactos con provenance exacta hace falta una ruta/ref autorizado por el owner que respete el contrato; hasta entonces, los cuatro paquetes del SHA `f0151cf…` son **NOT BUILT**. El preview anterior sirve como evidencia histórica de sus propios paquetes, no como entrega candidata.

La matriz `docs/release/manual-acceptance-0.3.0.md` no se modificó. Los casos gráficos/químicos que carecen de observación identificable permanecen `NOT TESTED`; los resultados manuales previos se conservan sin ampliarlos. No hubo instalación Windows/Flatpak, aceptación gráfica interactiva, sesión OpenGL, ni validación manual del propietario. Los smoke automatizados y el import OSTree temporal no cambian esos estados.

**Dictamen: NO-GO para tag/publicación beta.** Bloqueos: (1) faltan cuatro builds del SHA candidato exacto, (2) aceptación manual P0/P1 sigue pendiente, (3) CI verde #39782991976 tiene un shard con sólo 0.602 s de margen y el run previo del mismo código funcional tuvo timeouts 5/7, (4) las pruebas locales del validator PyInstaller completo no están disponibles. No se creó tag/Release, merge/PR, no se tocó `main`, `gh-pages` ni canales públicos.

### Verificación del registro

Tras actualizar documentación: `openspec validate final-stabilization-chemuson-0-3-0-beta-1 --strict` válido; `python -m compileall -q src tests tools packaging` y `git diff --check` pasan. `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` conserva el único F401 basal `math` sin usar en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó código ni se amplió el alcance para corregirlo. Se verificó que `docs/release/manual-acceptance-0.3.0.md` no tiene diff.

## Addendum — Build Preview RC exacto y descarga de candidata (2026-10-22)

### Referencia y ejecución de workflow

- Baseline inicial (registrada en `baseline.md`): checkout único `/home/ccachyavgp/Documentos/ChemUSON`, rama activa `stabilize/v0.3.0-beta.1`, HEAD `a38fa8045c858f8ff4b6002950f57e50292aea81`, igual a `origin/stabilize/v0.3.0-beta.1`, árbol limpio y sin worktrees adicionales. No se cambió el HEAD activo.
- Se publicó con push normal, sin `--force` ni commit nuevo, la referencia remota `release/v0.3.0-beta.1-rc-prep` apuntando al SHA exacto `a38fa8045c858f8ff4b6002950f57e50292aea81`. El ref se conserva; no se ha borrado ni movido. El HEAD/worktree activo siguió en `stabilize/v0.3.0-beta.1`.
- El push no inició un run de push; el commit existente contiene `[skip ci]`. El workflow declara `workflow_dispatch`, así que se lanzó desde **la misma rama RC**, sin cambiar source ni filtros: [Build Preview #39789295451](https://github.com/PJGV333/Chemuson/actions/runs/39789295451). Evento `workflow_dispatch`; `headBranch=release/v0.3.0-beta.1-rc-prep`; `headSha=a38fa8045c858f8ff4b6002950f57e50292aea81`. El job `Validate preview source` comprobó rama, versión canónica `0.3.0-beta.1` y SHA exacto; la validación local del mismo contrato también pasó.
- Conclusión final `success`, `complete=true`, `publication=false`; jobs de source validation, Windows portable + Inno, Linux AppImage Type 2, Linux Flatpak y reporte: todos `success`. Tiempos job: 5 s, 4m46s, 4m33s, 5m12s y 4 s respectivamente. El reporte está en [artifact 12430367969](https://github.com/PJGV333/Chemuson/actions/runs/39789295451/artifacts/12430367969). El event no creó Release, tag ni canal.

### Artifact groups y procedencia comprobada

Los cinco artifact groups (dos paquetes Windows producidos por un job) se descargaron desde Actions a `/tmp/chemuson-beta1-candidate-a38fa804/`; retención Actions configurada a 14 días. `sha256sum -c checksums.sha256` aprobó paquete y `preview-provenance.json` en los cuatro grupos. Para cada provenance se verificó versión, `build_type=preview`, `build_status=success`, SHA completo, rama RC y `publication=false`; tamaño y hash del paquete coinciden con provenance y checksum.

| Paquete | Artifact Actions | Archivo | Tamaño | SHA-256 del paquete | SHA-256 de `preview-provenance.json` | SHA-256 de `checksums.sha256` |
|---|---:|---|---:|---|---|---|
| Windows portable | [12429864673](https://github.com/PJGV333/Chemuson/actions/runs/39789295451/artifacts/12429864673) | `Chemuson-v0.3.0-beta.1-preview-a38fa804-windows-x86_64-portable.exe` | 133,249,529 | `816cce117e97dd4593e2980c9decaf70e7f10152cf019a3e74928d54f16418f4` | `6da7d2492dd83aea6bb66412c6668245b03a79c7a33331ee0d5b9b4df80a1d96` | `c95ea2f2221052dee5e90752a3a95fc3ff8a92f58103334ab4e062baaffaf5f3` |
| Windows Inno Setup | [12430132855](https://github.com/PJGV333/Chemuson/actions/runs/39789295451/artifacts/12430132855) | `Chemuson-v0.3.0-beta.1-preview-a38fa804-windows-x86_64-setup.exe` | 134,387,966 | `9550d40341401a1c1ef8c419b81790513281b84690cc647976c62610fb2b482d` | `1dd42f42cbe40e16c8b2a63ea81dfa4212390478a644c24faf4901e41318fc94` | `86be6bd5259f02aabdb7e6a7a7bad71b1a31a456bed6b44d39e935dc279d9d2b` |
| Linux AppImage Type 2 | [12431205949](https://github.com/PJGV333/Chemuson/actions/runs/39789295451/artifacts/12431205949) | `Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.AppImage` | 176,645,312 | `5559e28117d576b1b52249dc395aa42abd48b9cbc4bb2f84c483efdf8a643ea1` | `b119d8f2b46cfe176863e4866d36cbd3b084dca968c3b24a1e15bdbcfd400465` | `fa798f47f377eee449a90e071acef8cab03dc577c385135c9b2abd75ba676055` |
| Linux Flatpak | [12429843801](https://github.com/PJGV333/Chemuson/actions/runs/39789295451/artifacts/12429843801) | `Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.flatpak` | 110,084,016 | `8d83724b1795fc6d344bc2db162133b004a2498ff648b67454278ae1ee3d42d9` | `c5eba3b2b0307f5ec6ecbd49ba69b5044d73ffd0a76939af56704f0a46150fa7` | `e8e261de6db5d2f605b1002ff9515359da56aa000fc5767bd03a26d8a75977b3` |

El build report JSON declara los mismos cuatro SHA/branch/version/`publication=false`; tamaño 472 bytes, SHA-256 `2786c25c05cecc7e2e48acd2361afdbeb5aaf32f9e9188f1837b8b864420261d`. Se conserva junto con paquetes y logs locales en `/tmp/chemuson-beta1-candidate-a38fa804/`. Los URLs anteriores descargan artifacts ZIP de Actions que contienen paquete, checksum y provenance.

### Automatización sobre paquetes reales

**Windows portable.** `file` identifica PE32+ GUI x86-64. En el runner Windows, sobre el ejecutable que se copió al artifact portable, pasaron: Python ChemName acceptance pre-freeze; `validate_packaged_icons.py` (69 SVG visibles en tema claro/oscuro, 11 iconos/controles esenciales, DPR 2); `validate_packaged_rdkit_worker.py` (extensiones nativas `rdchem`, `rdMolDescriptors`, `rdDistGeom` bajo `_MEIPASS`, proceso padre sin RDKit/GUI, `CCO`, descriptors y 3D); y `validate_packaged_chemname.py` (9/9 templates, 6/6 moléculas comparadas con Python). El hash reportado por el smoke congelado ChemName `816cce117e97dd4593e2980c9decaf70e7f10152cf019a3e74928d54f16418f4` coincide con el portable descargado. No se ejecutó interactivamente en Windows en esta sesión.

**Windows Inno Setup.** Log registra compilador Inno Setup 6.7.1, compile exit 0 en 21.578 s y setup PE32 GUI Intel i386; el payload portable validado es x64. El artifact coincide con provenance/hash. No se instaló, actualizó ni desinstaló; aceptación del instalador queda NOT TESTED. PyInstaller emitió warnings de resolución de dependencias Qt opcionales (multimedia/PulseAudio, pcsclite, XCB, QML/3D) en el log; el workflow y las verificaciones requeridas acabaron en success. No se presenta esa salida como prueba manual de esas rutas opcionales.

**Linux AppImage.** El runner pasó `validate_appimage.py` completo sobre el paquete: ELF64 x86-64 Type 2 (`AI\x02`), AppDir/desktop/AppStream y CArchive/resources, 69 SVG, ChemName 9/9 + 6/6, RDKit/SMILES/descriptores/3D y headless startup de 12 s; ejecutable del bundle `--version` coincide con `0.3.0-beta.1`. El paquete descargado volvió a pasar localmente header/hash, extracción `--appimage-extract` sin FUSE, `AppRun --version`, y los validadores de icons, RDKit y ChemName contra `usr/bin/Chemuson`. Los reports locales confirmaron 69 iconos en cada tema, HiDPI DPR 2, `CCO` (logP -0.0014, TPSA 20.23, HBD/HBA 1), 9 templates y 6 nombres; sus JSON se conservan en el directorio de evidencia. `--appimage-updateinformation` es vacío y `.updateinfo`, `.update.json`, `.zsync` no existen. Pruebas de UI visual/install siguen NOT TESTED.

**Linux Flatpak.** El runner construyó bundle local y ejecutó en su `/app` sandbox-build smoke: 9 templates y 6 nombres ChemName; reporte sin errores. El artifact descargado se importó localmente a un OSTree temporal no instalado: ref `app/io.github.PJGV333.Chemuson/x86_64/preview-a38fa804`, commit `465e4509af9715349006ff4c0e8daa6c6e1fc72bf5500afea5603ae6d6cec5f0`. Checkout verificó App ID/launcher, AppStream versión, y 9/9 templates `.mol`. No se instaló Flatpak en perfil de usuario/sistema, no se añadió remoto ni se inició su GUI.

### Evidencia local, manual y alcance

- Captura completa del log Actions: `/tmp/chemuson-beta1-candidate-a38fa804/actions-full.log`, SHA-256 `1c891d0f8fb4cc1617972bc13b847488f2c7971b8cccc4cd50789ea85ac301b7`. Descargas, manifests, JSON de las tres comprobaciones AppImage, extracción y checkout OSTree permanecen en el mismo directorio externo; manifiesto `EVIDENCE-MANIFEST.sha256` cubre y verifica 44 resultados/paquetes (SHA-256 del manifiesto: `d8672314b6eed69b28197521af6f92919af56a14f7ef582599163e4adb610348`). No se añadieron binarios al repo.
- La matriz `docs/release/manual-acceptance-0.3.0.md` no se cambió. Windows portable/startup/dibujo/guardado/ChemName/RDKit/icons visuales, instalación/upgrade/uninstall del setup, AppImage GUI manual y Flatpak instalado/GUI siguen NOT TESTED. No se marca PASS manual ni de capturas inexistentes.
- El repo no cambió código/workflow ni dependencias; catálogo de arquitectura intacto. No se hizo merge, PR, tag, GitHub Release, publicación ni cambio a main/canales. La referencia remota `release/v0.3.0-beta.1-rc-prep` debe mantenerse en el SHA exacto hasta que el owner complete la aceptación manual.
- Guía concisa de verificación e instalación para CachyOS/Windows: `docs/release/PREVIEW_BUILDS.md`. Dictamen: **candidata trazable lista para aceptación manual; no publicar aún**. Los P0/P1 manuales continúan bloqueando cualquier beta pública.
- Se verificó la procedencia exacta del build y los hashes; no se lanzó una segunda compilación independiente, por lo que no se afirma reproducibilidad binaria bit-a-bit entre ejecuciones.
- Cierre del checkout: HEAD local sigue exactamente `a38fa804…`, la referencia RC remota coincide y no se creó commit. Para conservar las actualizaciones pedidas de OpenSpec/release hay seis archivos documentales con cambios unstaged; por tanto `git status --short` no está vacío. Sin commit/reset/hide no se pueden conservar esas ediciones y tener a la vez árbol Git limpio. No hay binarios ni temporales generados dentro del repo.

## Addendum — integración, auditoría de seguridad y cierre público (2026-10-23)

### Integración autorizada

El propietario pidió expresamente cerrar la integración en `release/v0.3.0-beta.1-prep`. Tras `git fetch origin`, el ancestro común fue `61f7647b683db3eca19ea06e5ec82ab2897d640a`; stabilization estaba 29 commits adelante y 0 detrás. Árbol limpio, sin archivos sin merge ni operaciones Git incompletas. `git merge --ff-only stabilize/v0.3.0-beta.1` llevó release-prep al SHA `d808070fcf8e9efb1caf2e8a3fb6f663653b69a1`; los árboles eran idénticos. Push normal a origin completado; `origin/release/v0.3.0-beta.1-prep` coincide. `main` sigue en `f0fde603371255902bf0e630ca1ce032b8f16ad3`; no hubo force-push, tag, Release ni cambio de canal.

La primera ejecución automática posterior a ese push, Build Preview [#39893346523](https://github.com/PJGV333/Chemuson/actions/runs/39893346523), validó el source y luego se canceló deliberadamente al preparar la revisión pública del cuerpo de release. No produjo ni se utilizará ningún paquete. El preview candidato debe ser el nuevo run automático después del commit documental, no #39893346523 ni los artifacts históricos de `7d935669…`.

### Matriz compacta de presencia

| Área | Presencia en código y commits representativos | Evidencia / alcance de aceptación |
|---|---|---|
| Editor 2D: enlaces/anillos, etiquetas y continuidad | `src/chemuson/gui/canvas/{canvas_bond_geometry.py,canvas_structure.py,canvas_text.py}`; geometría inicial `32fc76f5`; amidas/ureas `3be1d700`, `582fd0e2`; aminas primarias/secundarias `a7fd1206`, `a5336866`. | `tests/test_amine_label_orientation.py`, pruebas de geometría/etiquetas y evidencia clara/oscura en `openspec/changes/fix-primary-amine-label-orientation-beta1/evidence/`. El propietario aceptó visualmente NH₂/H₂N horizontales y NH/HN en sector abierto, preservando amidas/ureas; esa aceptación es de la corrección/capturas, no de paquetes finales. |
| Clean2D, edición y persistencia | `src/chemuson/clean2d/`, `src/chemuson/gui/commands/`, `src/chemuson/chemio/persistence.py`; campañas `9d28c457`, `23b70462`, `946de072`; round-trip y exportación MOL `34711976`, `6df2ac01`. | Corpus de Clean2D y pruebas `tests/test_mol_interoperability_p0.py`, `tests/test_persistence_contract.py`, `tests/test_geometry3d_persistence.py`; CI #39782991976 del SHA anterior registró cobertura automatizada. La aceptación manual amplia de guardar/reabrir, importación y recuperación continúa `NOT TESTED`; `DATA-11` sólo tiene un PASS limitado a una molécula. |
| Visor 3D, conformeros, worker, edición y Live MM | `src/chemuson/gui/viewport3d/`, `src/chemuson/geometry3d/`, `src/chemuson/chemio/rdkit_safe.py`; geometría/editor `d8f8a2d2`, `f78c58fc`, Live MM `6e53b408`, profundidad `ec486626`, worker congelado `6a6342f9`. | El propietario aceptó visualmente el AppImage diagnóstico en CachyOS (2D y 3D, incluida cafeína). Build Preview #39881978237 en SHA `7d935669…` pasó 2D→3D congelado en Windows/AppImage. `geometry_flow_smoke` de fuente también pasó en `d808070f…`; el preview exacto de ese SHA fue cancelado y no sirve como aceptación. Windows/Flatpak/Live MM manual continúan pendientes. |
| UI, menús, paneles y herramientas | `src/chemuson/gui/main_window_ui_builder.py`, `side_panel`, `viewport3d`, `energy_diagrams.py`; contrato de Ver/3D `409581dc`; rótulos de energía `61f7647b`; Assistant en controladores y `src/chemuson/molecular_assistant/`. | `tests/test_beta_compact_shell_and_tabs.py` y `tests/test_side_panel.py` cubren el menú Ver/visor 3D y los siete docks del panel lateral más el octavo dock 3D independiente; `tests/test_energy_diagrams.py` y los tests del Assistant cubren contratos acotados. ChemName está en `src/chemuson/chemname/`; tests source incluyen `tests/test_chemname_acceptance.py`. Assistant es experimental y requiere revisión; no se presenta cobertura general IUPAC como aceptación. |
| Distribución y política de builds | `chemuson.spec`, `packaging/release/`, `packaging/linux/`, `.github/workflows/{build-preview.yml,release.yml,test.yml}`; aislamiento AppImage/worker `db48d715`, `6a6342f9`; pruebas locales/salida externa `201a81cb`. | Preview #39881978237 compiló los cuatro formatos desde su SHA exacto con checksums/provenance y smoke 2D→3D. **No es el SHA final**: las correcciones posteriores de aminas aún requieren preview exacto. `CHEMUSON_OUTPUT_DIR` está implementado. Sin paquete nuevo final: los cuatro formatos actuales quedan `NOT BUILT` hasta el siguiente run. |

La matriz demuestra presencia de código y evidencia asociada; no transforma pruebas source, screenshots ni previews de SHA anteriores en aceptación manual del candidato exacto.

### Auditoría focalizada de integridad y seguridad

- **Archivos/datos:** `FileController.open_file_path` abre en pestaña nueva y la descarta ante parseo fallido; el documento previo queda separado. `.cmsn` usa JSON y escritura temporal seguida de `os.replace`; importadores MOL/CML no evalúan contenido como código. Se revisaron persistencia, controladores y rutas de entrada/salida. No se confirmó traversal ni ejecución de código por archivo. No hay límites generales de tamaño para todos los formatos, por lo que archivos enormes pueden consumir recursos; no se demostró pérdida de datos ni se cambió el formato.
- **Procesos/workers:** los usos de `subprocess` revisados usan vectores de argumentos, límites temporales y no `shell=True`; el worker RDKit usa `sys.executable`/modo interno y JSON temporal, y Open Babel usa `shutil.which` más argumentos separados. No se encontró `eval`, `exec`, `pickle.loads` ni carga YAML insegura en las rutas inspeccionadas. Esto es revisión estática focalizada, no pentest.
- **Assistant/secretos:** `src/chemuson/molecular_assistant/provider.py` valida HTTP(S), prohíbe credenciales/query en URL, no sigue redirects, limita respuesta/timeout y oculta errores/secretos. Claves se envían sólo por HTTPS o HTTP loopback y `platform/settings.py` no las persiste (elimina la clave legacy). El endpoint HTTP remoto sin clave sigue permitido intencionalmente por OpenSpec; para datos privados se debe usar HTTPS. Los prompts/estructuras enviados al Assistant remoto quedan sujetos al proveedor. La consulta externa de identidad requiere opt-in y se desactiva por defecto.
- **Updater:** `src/chemuson/update/` exige URL HTTPS, tiene límites de descarga, `require_sha256=True` por defecto, control semver/canal y rollback; firma criptográfica sigue opcional, no se afirma que todos los futuros artifacts estén firmados. No se confirmó bypass explotable.
- **Dependencias:** `pip-audit -r requirements.txt -r requirements-dev.txt` resolvió 13 dependencias, Pillow `13.0.0`, **sin vulnerabilidades conocidas** al 2026-10-23. El entorno local separado contiene Pillow `12.3.0`, para el que el escaneo informa CVE-2026-96637, CVE-2026-96643, CVE-2026-96644, CVE-2026-96778 y GHSA-78g4-xqwv-q3g8. Ningún módulo bajo `src/` importa PIL y `chemuson.spec` no lo declara como hidden import; el entorno local no representa el runtime resuelto del build. Falta confirmar el inventario CArchive del nuevo paquete. Dependabot alerts están deshabilitadas y no hay análisis Code Scanning disponible (`no analysis found`); estos controles quedan **NOT TESTED**, no PASS.
- **Hallazgos P0:** no se confirmó en el código revisado corrupción/pérdida provocada, exposición de claves, ejecución de datos de archivos como código o vulnerabilidad grave alcanzable. La auditoría no es una prueba de penetración. Los casos manuales P0 de apertura/cierre, guardar/reabrir, importación e identidad estereoquímica siguen pendientes en paquetes exactos.

### Bloqueos históricos P0/P1 y evidencia actual

| Hallazgo histórico | Estado frente al código y última evidencia | Decisión actual |
|---|---|---|
| P0 persistencia/importación/estereoquímica/cierre | Tests source focales pasan; el propietario sólo confirmó un round-trip MOL estereoquímico limitado y el smoke visual 3D diagnóstico. | Las filas restantes (datos corruptos, identidad, cierre activo, save/reopen) son `NOT TESTED`; no hay P0 confirmado, pero aceptación esencial pendiente bloquea tag. |
| P1 iconos ausentes en paquete | Fallo histórico #37826597134; los gates empaquetados posteriores comprobaron SVG visibles en sus propios SHA. | Repetir/verificar el smoke y que el owner inspeccione los paquetes exactos Windows y Linux; no convertir el fallo histórico en PASS por inspección source. |
| P1 worker RDKit/viewport 3D congelado vacío | Fix `6a6342f9`; aceptación visual del propietario en diagnóstico y smoke congelado en preview #39881978237. | Code fix confirmado; nueva verificación automática de paquete para el SHA final más retest GUI compatible siguen pendientes. |
| P1 templates ChemName ausentes | Falla histórica de CArchive corregida; posteriores previews ejecutaron 9/9 templates y 6/6 controles. | Repetir gate en el nuevo build; ChemName GUI del candidato aún no está aceptado manualmente. |

### Protección de publicación y decisión

Consultas GitHub con sesión autenticada: `rulesets` devuelve `[]`; la rama prep no está protegida; no existe protección `v*` observable (legacy endpoint no disponible y ningún ruleset); los únicos entornos configurados son `github-pages`, sin reviewers. Los entornos dinámicos `beta` y `flatpak-beta` no están configurados. `default_workflow_permissions=read`; sólo los jobs oficiales `release` y `publish_flatpak_remote` solicitan `contents: write`, pero `allowed_actions=all` y `sha_pinning_required=false`. El workflow de release sigue tag-only (`v*`), gate acotado, build desde SHA validado, `prerelease=true` en canal beta y `overwrite_files=false`; ahora incluye el cuerpo curado de `docs/release/0.3.0-beta.1-public.md` mediante `body_path`.

**NO-GO para tag/publicación:** (1) Build Preview exacto y paquete CArchive no disponibles todavía; (2) aceptación humana P0/P1 sobre los cuatro artifacts exactos sigue `NOT TESTED`; (3) faltan protección de tag `v*` y entornos beta protegidos. Se preparó una lista corta para el propietario en `docs/release/manual-acceptance-0.3.0.md`. No crear tag, release ni desplegar Flatpak hasta resolver estos gates. `main` permanece intacta.

### Política de pruebas y resultado local

`AGENTS.md` exige históricamente suite completa en baseline/final; `docs/testing/AGENT_TESTING_POLICY.md` posterior declara que pytest masivo es extraordinario. Se registró explícitamente la contradicción y una propuesta para reconciliarla, sin modificar `AGENTS.md`, borrar pruebas o reactivar la suite. Baseline fuente `d808070f…`: compileall PASS, colecta 2.624, arquitectura 283/283, conjunto focal 220/220, smoke fuente real CCO→RDKit 3D PASS, OpenSpec strict PASS. Ruff focal conserva sólo el F401 basal `math`. La ejecución completa no se lanzó por la política específica y solicitud del propietario.

El resultado final del Build Preview que se ejecute después del commit documental se reporta con el SHA, URLs, artifact IDs y hashes en la respuesta al propietario. No se añadirá un commit de resultados después de ese run porque alteraría el SHA fuente que los paquetes deben acreditar; los manifests y el reporte de Actions son la fuente verificable de esa ejecución.
