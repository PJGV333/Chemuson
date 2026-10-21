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