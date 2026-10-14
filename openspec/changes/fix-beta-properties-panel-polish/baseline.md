# Baseline — microcampaña final del panel

Captura antes de tocar código o documentación del producto. Rama creada desde el HEAD solicitado; sólo los artefactos OpenSpec de esta campaña se escriben después de esta captura.

- Base indicada/verificada: `release/v0.3.0-beta.1-prep` en `8f20e675d8be8632c5d750e63699830d6d1223bb`; `git status --short` vacío.
- Rama de trabajo: `fix/beta-properties-panel-polish`, HEAD inicial `8f20e675d8be8632c5d750e63699830d6d1223bb`.
- `python -m compileall src tests tools packaging` (timeout externo 300 s): exit 0, `ELAPSED=0.20 s`.
- `pytest --collect-only -q` (timeout externo 600 s): exit 0, `2147 tests collected in 0.87s`, tiempo total `1.27 s`.
- `pytest -q` (timeout externo 600 s): exit 1, `1 failed, 2126 passed, 20 skipped in 595.74s`, tiempo total `596.81 s`. Único fallo: `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`; el orden observado fue `[('worker.finished', 1)]` y faltaron `thread.finished`/`job_finished`. Es un fallo previo al cambio en un área fuera de alcance; no modificar CompChem. No se volverá a ejecutar el monolítico.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` (timeout externo 300 s): exit 1 sólo por el baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no tocar.
- El baseline completo tardó 9m56.81s y estuvo bajo el timeout externo de 10 minutos solicitado. Las validaciones posteriores serán shards con timeout externo ≤600 s, preferentemente 120 s para Qt.

## Contratos revisados

- Cambio OpenSpec activo: `prepare-chemuson-0-3-0-beta-1` (proposal/design/tasks/specs leídos); specs aplicables `ui-side-panel-status-bar`, `ui-polish` y delta de onboarding/aceptación.
- El SidePanel actual fija 340 px, se añade después del lienzo en un `QHBoxLayout`, inicializa visibilidad de `ui/side_panel/visible` y conserva correctamente `active_tab`; el menú Ver ya tiene siete acciones directas.
- `ChemicalPropertiesDock` conserva `_rows` fuente y copia TSV, pero renderiza cada celda con texto plano.
- El onboarding actual usa SidePanel como tercer target aunque el panel puede estar oculto.
