# Validación

## Baseline

- Rama de trabajo: `fix/beta-save-shortcut-menu-access`; HEAD inicial `0765dcba7abf999e5f3c56efa5634d55bc2d2348`, árbol limpio antes de esta corrección.
- `python -m compileall src tests tools packaging`: PASS en baseline y tras los cambios.
- Baseline `pytest --collect-only -q`: 2169 IDs; el manifiesto actualizado cubre 2174 IDs.
- Baseline monolítico `pytest -q`, timeout externo 480 s: venció el límite (exit 124); antes emitió un fallo transitorio en `test_compchem_controller_generates_async_with_fake_backend`, que pasó al ejecutarse aislado. No se modificó esa área.
- Ruff baseline con `--select F401,F811,F821,E722,E741`: el único hallazgo preexistente es `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se alteró fuera del alcance.

## Verificación de esta corrección

- Suite completa mediante ocho shards acotados, cada uno con timeout de 300 s. `tools/ci_pytest_shards.py verify` confirmó **2174/2174 IDs únicos**, los ocho shards completos, exit 0, sin timeout ni fallos: **2154 passed, 20 skipped**.
  - Shard 0: 272 passed.
  - Shard 1: 269 passed, 3 skipped (repetido tras los ajustes finales de imports del test).
  - Shard 2: 270 passed, 2 skipped.
  - Shard 3: 267 passed, 5 skipped (repetido tras los ajustes finales de imports del test).
  - Shard 4: 272 passed (repetido tras retirar una asignación no usada del test).
  - Shard 5: 269 passed, 3 skipped.
  - Shard 6: 268 passed, 3 skipped.
  - Shard 7: 267 passed, 4 skipped.
- Pruebas focalizadas de UI y preferencias: PASS; arquitectura: **280 passed**; contrato del workflow de pruebas: **4 passed**.
- `python -m compileall -q src tests tools packaging`: PASS.
- Ruff para los módulos Python cambiados, reglas `F401,F811,F821,E722,E741`, e imports/variables locales de pruebas cambiadas, reglas `I,F841`: PASS. El Ruff de todo el árbol conserva únicamente el F401 baseline indicado arriba.
- `git diff --check`: PASS.
- `openspec validate fix-beta-tab-close-and-toolbar-visibility --strict --no-interactive`: PASS.
- El Ctrl+S probado sigue llegando a la acción y ruta de guardado existentes; no se cambió la serialización ni el flujo CMSN.

## Aceptación y publicación

- La aceptación visual interactiva en las plataformas del propietario queda pendiente; pruebas Qt/render y automatización no se presentan como sustituto de esa aceptación manual.
- Commit de implementación: `a014650cc8a0c215a7276159de7a41374433ed1c` (`fix(gui): show tab close controls and main toolbar`), publicado por push normal únicamente a `fix/beta-save-shortcut-menu-access`.
- GitHub Actions [#38779862340](https://github.com/PJGV333/Chemuson/actions/runs/38779862340): pytest-plan, Windows-smoke y Flatpak-smoke pasaron; la validación de pytest quedó roja. El shard 1 falló dos veces (intento 1 y reejecución de fallidos) únicamente en el caso baseline `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`: a los 3 s el test había registrado `worker.finished`, pero aún no `thread.finished`/`job_finished`. El mismo nodo ya había fallado transitoriamente durante el baseline previo a cambios y pasó aislado; los shards locales finales pasaron. No se atribuye a este cambio de GUI.
- La CI verde queda bloqueada por este fallo preexistente e intermitente. No se modificó dicho test ni otro subsistema fuera del alcance autorizado; queda pendiente decisión sobre autorizar un endurecimiento test-only o aceptar/documentar el baseline existente.
- No se integra a beta/main ni se crean tags/releases.
