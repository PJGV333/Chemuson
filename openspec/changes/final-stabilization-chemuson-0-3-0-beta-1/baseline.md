# Baseline — estabilización final beta.1

Fecha: 2026-10-21 UTC. Rama de campaña local: `stabilize/v0.3.0-beta.1` (mismo directorio; sin worktree). Rama de origen: `release/v0.3.0-beta.1-prep`. SHA esperado/inicial: `61f7647b683db3eca19ea06e5ec82ab2897d640a`.

## Git / arquitectura / publicaciones

- Antes del cambio, `git status --short --branch` en la rama prep mostró `## release/v0.3.0-beta.1-prep` sin modificaciones; se confirmó HEAD esperado.
- Se creó únicamente la rama local `stabilize/v0.3.0-beta.1` desde ese SHA. HEAD de baseline sigue siendo `61f7647b683db3eca19ea06e5ec82ab2897d640a`; el árbol seguía limpio tras crear la rama.
- Se leyeron `AGENTS.md`, `architecture/modules.yml`, el OpenSpec `prepare-chemuson-0-3-0-beta-1` y el OpenSpec de CI reproducible, workflows `test.yml`, guías `docs/release/*` y matriz manual existente. No hay cambio de paquete/import/módulo previsto; no se requiere tocar el catálogo de arquitectura.
- No se creó tag, GitHub Release, PR, preview adicional ni push. No se modificó `main`, la rama prep, `gh-pages` ni canales públicos.

## Baseline local antes de editar

- `python -m compileall src tests tools packaging`: **exit 0**.
- `pytest --collect-only -q` con `/usr/bin/python` (Python 3.14.8 de sistema): **exit 2**, dos errores de colección por falta de `rdkit` en el intérprete (`tests/test_mol_interoperability_p0.py`, `tests/test_smiles_stereo_import.py`). Esto es una limitación de selección del intérprete, no una falla de código; no se usa para representar la colección del proyecto.
- `chem/bin/python -m pytest --collect-only -q`: **2.597 tests collected in 0.82s**, exit 0. Esta es la colección con el entorno local del proyecto (Python 3.14.8, pytest 9.1.1, PyQt 6.11.0/Qt 6.11.0, RDKit instalado).
- `timeout --signal=TERM --kill-after=10s 8m chem/bin/python -m pytest -q`: **timeout exit 124** cerca del 30% de progreso, sin resumen final. No se interpreta como PASS ni como una de las fallas CI; se conserva `/tmp/chemuson-final-stabilization-baseline/pytest-project-python.txt`. Se usarán los ocho shards bounded existentes, no se repetirá como monolítico.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: **exit 1** sólo por el baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; sin cambios.
- Salidas exactas de captura local en `/tmp/chemuson-final-stabilization-baseline/` (status, compileall, collect en ambos intérpretes, pytest bounded y Ruff). El log raw de colección del entorno de proyecto contiene la lista completa de 2.597 IDs.

## Run Actions que motiva esta campaña

- Run [#39679607676](https://github.com/PJGV333/Chemuson/actions/runs/39679607676), SHA exacto `61f7647b683db3eca19ea06e5ec82ab2897d640a`, rama prep; Python 3.11.17 / Ubuntu 24.04 / `QT_QPA_PLATFORM=offscreen`; conclusión `failure`.
- Plan: 2.597 IDs, 8 shards, distribución `[301,243,327,362,449,267,306,342]`, hash de manifiesto esperado `31c2db63a18ff6dae6877484ad95057d4a93bc59a28ccf3825940e599579dd95`.
- JUnit verificó las asignaciones exactas de los ocho shards; totals: **2.572 passed, 22 skipped, 3 failed, 0 errors**. Ningún shard tuvo timeout/señal de proceso ni IDs omitidos. Shards 0 y 6 fallaron, summary falló correctamente; `windows-smoke` y `flatpak-smoke` pasaron.

| Shard | IDs/JUnit | Passed | Skipped | Failed | Errors | Tiempo | Estado |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 301/301 | 298 | 2 | 1 | 0 | 178.216 s | failure |
| 1 | 243/243 | 241 | 2 | 0 | 0 | 122.944 s | success |
| 2 | 327/327 | 324 | 3 | 0 | 0 | 139.013 s | success |
| 3 | 362/362 | 360 | 2 | 0 | 0 | 173.152 s | success |
| 4 | 449/449 | 448 | 1 | 0 | 0 | 86.899 s | success |
| 5 | 267/267 | 267 | 0 | 0 | 0 | 219.351 s | success |
| 6 | 306/306 | 304 | 0 | 2 | 0 | 101.179 s | failure |
| 7 | 342/342 | 330 | 12 | 0 | 0 | 84.487 s | success |

- Artefactos Actions descargados y conservados temporalmente en `/tmp/chemuson-actions-39679607676/artifacts/`: plan artifact `12383167003`; shard artifacts 0 `12382874613`, 1 `12381989909`, 2 `12382074953`, 3 `12384190464`, 4 `12383815513`, 5 `12383462150`, 6 `12384410118`, 7 `12383437276`. Contienen JSON de proceso, JUnit, logs de pytest y timings JSONL.
- Falla A: `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`. El log registra sólo `worker.finished`; faltan `thread.finished` y `job_finished`, con `active_jobs=(1,)`, `pending_results=(1,)` después del límite fijo existente de 3 s. No es un timeout/abort del shard. Tres ejecuciones aisladas y el grupo local `tests/test_command_palette.py` + este node pasan en el host Python 3.14/PyQt 6.11; la discrepancia de runtime no descarta el fallo reproducido en Actions.
- Falla B: `tests/test_gui_async_worker_shutdown.py::test_close_waits_for_all_descendant_worker_families_and_suppresses_late_ui` lanza `AttributeError: 'HeldCompChemWorker' object has no attribute 'progress_ready'` en `CompChem3DController.start_job`, línea 183. El test falla antes de sus aserciones de cierre/espera y no valida ese comportamiento.
- Falla C: `tests/test_p0_qt_shutdown_subprocess.py::test_existing_async_worker_shutdown_contract_in_isolated_pytest_process` ejecuta el mismo archivo en un proceso hijo y propaga la falla B (`.F`, 1 failed, 1 passed). No es una tercera causa.
- El warning `QOpenGLWidget is not supported on this platform` es coherente con offscreen y no explica los failures. En este run no hubo SIGSEGV de pytest ni crash de Qt reportado por el supervisor.
- Salidas completas de run: artifacts enlazados por IDs anteriores; resumen agregado y confirmación de jobs descargados en `/tmp/chemuson-actions-39679607676/`.

## Estados de skips

Los 22 skips pertenecen al resultado automático preexistente del run. Antes de declarar aceptación CI, se identificarán los 22 node IDs y sus razones desde JUnit; no se aceptará ningún skip nuevo ni se actualizará el manifiesto para esconder fallos. Su detalle se añadirá a `validation.md`.

## Baseline de la candidata RC — 2026-10-22 (antes del push de referencia)

- Directorio único `/home/ccachyavgp/Documentos/ChemUSON`; rama local activa `stabilize/v0.3.0-beta.1`; HEAD `a38fa8045c858f8ff4b6002950f57e50292aea81`; `origin/stabilize/v0.3.0-beta.1` coincidía exactamente. `git status --short` no produjo salida (árbol limpio). `git worktree list` confirmó sólo el checkout principal. `git ls-remote` confirmó que `release/v0.3.0-beta.1-rc-prep` aún no existía antes del push.
- Diff de código desde el SHA del último CI verde (`f0151cf09a7200acc53cb6b371105e79fd3cac1d`) hasta `a38fa804…`: vacío en `src`, `tests`, `tools`, `packaging` y `chemuson.spec`; el cambio fue documental. `python packaging/release/preview_build.py verify --version 0.3.0-beta.1 --git-sha a38fa8045c858f8ff4b6002950f57e50292aea81 --source-branch release/v0.3.0-beta.1-rc-prep` devolvió `Verified preview source a38fa8045c858f8ff4b6002950f57e50292aea81` (exit 0); el validator acepta el nombre de rama RC dentro de `release/**-prep`.
- `git status --short`: salida vacía, exit 0. Captura exacta en `/tmp/chemuson-beta1-candidate-a38fa804-baseline/git-status-short.txt` (SHA-256 del archivo vacío: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`).
- `python -m compileall src tests tools packaging`: exit 0. Salida exacta preservada en `/tmp/chemuson-beta1-candidate-a38fa804-baseline/compileall.txt` (SHA-256 `6ec2f95b06698331cfb9508c8ab1a324349b4061e2bf81f891aa8e531fa3859e`).
- `timeout --signal=TERM --kill-after=5s 120s chem/bin/python -m pytest --collect-only -q`: exit 0; `2598 tests collected in 0.68s`. Salida completa de IDs en `/tmp/chemuson-beta1-candidate-a38fa804-baseline/pytest-collect.txt` (SHA-256 `2065e7d02bc6c0f8c8c5e6c78a3a50adbdc900ffad0fd42f2eb007bffadb09b0`).
- `pytest -q` completo no se repitió: el baseline monolítico exacto ya está registrado arriba, expiró a los 8 min sin resumen, y la guía bounded prohíbe repetirlo. El último CI completo disponible es #39782991976 sobre `f0151cf…`; respecto a `a38fa804…` sólo variaron documentos, sin código/tests/packaging. No se usó un shard local ni se lanzaron los ocho shards.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 únicamente por el F401 basal `math` no usado en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. Captura exacta `/tmp/chemuson-beta1-candidate-a38fa804-baseline/ruff.txt` (SHA-256 `1f8ffaf9f6a995877e4fa56605f8ad9d07f0499ef11fb3787e3c6619ff22cf19`); no se corrigió fuera de alcance.
- Capturas completas del baseline y preflight viven en `/tmp/chemuson-beta1-candidate-a38fa804-baseline/`; no se copiaron binarios ni datos de build al repositorio.
