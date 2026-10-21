# Validación — disposición visual N–H

## Gates del cambio

- Tests dirigidos: `QT_QPA_PLATFORM=offscreen .venv/bin/pytest -q tests/test_amide_urea_2d_continuity.py tests/test_functional_group_hydrogen_adjustment.py` — **32 passed**. Incluyen eventos Qt reales para urea, acetanilida y aminas, geometría, anclaje del N, temas, zoom/DPR, undo/redo y CMSN.
- Matriz visual/funcional focal ampliada — **444 passed** (capturas y métricas en `evidence/`).
- `python -m compileall src tests tools packaging` — PASS.
- Ruff focal (`F401,F811,F821,E722,E741`) en los cuatro archivos Python modificados — PASS.
- `openspec validate fix-amide-urea-nh-label-layout-beta1 --strict` — PASS.
- `git diff --check` y revisión de whitespace de Markdown no rastreado — PASS.
- Capturas comparativas Qt offscreen claras/oscuras para urea, acetanilida y aminas. No equivalen a aceptación visual en escritorio nativo.

## Suite completa final

Los 2.597 node IDs del manifiesto se ejecutaron mediante ocho shards LPT, cada uno con límite externo de 300 s, en serie para evitar contención de CPU. JUnit y timings cubren exactamente la asignación en los ocho shards. Resultado agregado: **2.572 passed, 22 skipped, 3 failed**. El verifier confirmó cobertura completa, pero correctamente rechazó el estado por fallos en shards 0 y 6. No hubo timeout en la ejecución final.

Dos causas preexistentes explican los tres node IDs fallidos:

1. `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` — falla solo en el grupo ordenado de Qt, dejando únicamente `worker.finished` y sin completar la entrega diferida. Reproducido con los primeros 125 node IDs ordenados del shard 0 en el worktree base inmutable `409581dc685fd166428bde47d8355633a8eb05da`: **123 passed, 1 skipped, 1 failed**. El nodo pasa al ejecutarse aislado. No se modificó ese test ni su controller.
2. `tests/test_gui_async_worker_shutdown.py::test_close_waits_for_all_descendant_worker_families_and_suppresses_late_ui` — su fake `HeldCompChemWorker` no declara `progress_ready`, que el controller conecta. Reproducido directamente en el mismo commit base. `tests/test_p0_qt_shutdown_subprocess.py::test_existing_async_worker_shutdown_contract_in_isolated_pytest_process` vuelve a ejecutar ese test y reproduce el mismo error; por eso una sola incompatibilidad origina dos node IDs fallidos. Ninguno de esos archivos fue modificado.

La primera ejecución concurrente de shards tuvo contención: shard 5 alcanzó el límite y una prueba UI falló de forma transitoria. Ambos se reejecutaron en serie dentro del límite de 300 s: shard 5 terminó en 247,669 s y la prueba UI pasó. Esos resultados transitorios no forman parte del agregado final.

El baseline monolítico previo agotó el timeout de 600 s alrededor del 30%; por eso la suite final se completó y auditó con shards. Los fallos anteriores no se ocultaron ni se arreglaron fuera del alcance. El propietario autorizó explícitamente el cierre con commit y push pese a esos fallos basales, que se conservan documentados; no se hizo excepción en el verifier ni se tocaron esos tests. No se autorizaron merge, PR, tag ni Release.
