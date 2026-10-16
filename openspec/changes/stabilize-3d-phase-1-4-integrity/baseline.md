# Baseline — stabilize-3d-phase-1-4-integrity

Fecha de captura: 2026-10-16. Capturado antes de modificar código o documentación del cambio.

- Rama: `feature/chemuson-3d`; upstream `origin/feature/chemuson-3d`.
- HEAD: `e53e7a23362359e66730b2af2645f7c51b6c6ad6` (`docs(openspec): record derived hydrogen continuity`), igual al upstream.
- Worktree inicial: `git status --short` vacío.
- Entorno: `.venv` del repo; CPython 3.13.14, pytest 9.1.1, Ruff 0.14.14; RDKit instalado.

## Salidas de baseline de esta continuación

- `./.venv/bin/python -m compileall src tests tools packaging`: exit 0; salida guardada en `/home/ccachyavgp/.hermes/cache/scratch/chemuson-campaign-baseline-e53e7a2/compileall.log`.
- `./.venv/bin/pytest --collect-only -q`: exit 0; 2150 tests collected in 0.66s; lista completa guardada en `/home/ccachyavgp/.hermes/cache/scratch/chemuson-campaign-baseline-e53e7a2/collect.log`.
- `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por un hallazgo previo fuera del alcance: `F401` `math` sin uso en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; salida exacta en `/home/ccachyavgp/.hermes/cache/scratch/chemuson-campaign-baseline-e53e7a2/ruff.log`.
- Pruebas focalizadas pre-cambio: `./.venv/bin/pytest -q tests/test_geometry3d_service.py tests/test_rdkit_live_stream.py tests/test_compchem3d_dock.py tests/test_geometry3d_editing.py tests/test_geometry3d_persistence.py tests/test_viewport3d_widget.py tests/test_viewport3d_camera.py tests/test_viewport3d_picking.py tests/test_gui_async_worker_shutdown.py`: **80 passed, 1 skipped in 21.73s**.
- La suite completa `pytest -q` no se repitió: los informes previos de esta rama conservan la salida/traceback del SIGSEGV baseline en `concurrent.futures.thread._worker`, y el encargo explícitamente pide no repetirla sin necesidad diagnóstica concreta. No se toma la omisión como un pase.

## Evidencia previa conservada (no reescrita)

- Baselines y logs históricos de `add-interactive-3d-molecular-viewport` y `add-editable-3d-molecular-editor` permanecen sin cambios.
- `validation.md` de ambos cambios mantiene sus resultados previos: tests, demo de butano, límites de autorización y crashes documentados.
- El core de la demo Qt/SIP se analizó de manera no interactiva desde `coredumpctl` y se copió a scratch; GDB corroboró `sip_api_get_address(w=0x190)` → `cleanup_qobject` → `sip_api_visit_wrappers` → `release_QApplication` en el hilo principal al finalizar el script. GDB advirtió que el core podría no coincidir exactamente con el ejecutable suministrado (se pasó el symlink `.venv/bin/python`); los módulos PyQt6/Qt de la traza coinciden con el entorno `.venv`. No se obtuvo causa raíz. Core extraído: `/home/ccachyavgp/.hermes/cache/scratch/chemuson-phase4-core-1490286`.

## Alcance del baseline

La colección, compilación y pruebas focalizadas sí fueron ejecutadas en HEAD `e53e7a2`. No se afirma que el SIGSEGV Qt/SIP o el cierre de display se haya reproducido en esta captura; la comparación real en display permanece bloqueada y reservada para prueba manual.