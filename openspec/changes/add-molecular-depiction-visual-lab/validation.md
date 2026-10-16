# Validación — Molecular Depiction Visual Lab

## Alcance ejecutado

Rama `feat/molecular-depiction-visual-lab`, creada desde `ee5294800c5779f7451e625d8089025e8e6632ab`; OpenSpec inicial quedó separado en `1af8f7e48e36472f62b688434f20b908baf82522`. No se cambiaron Clean2D, ChemName, ChemIO, `.cmsn`, dependencias, preferencias de la aplicación, `architecture/modules.yml`, ramas protegidas ni canales de publicación. El código de dibujo nuevo es opt-in; todos los overrides de `DrawingStyle` tienen `None` por defecto y `Current` mapea exactamente a `CHEMDOODLE_LIKE`.

## Verificación focal

- `pytest -q tests/test_depiction_lab_parameters.py tests/test_depiction_lab_renderer.py tests/test_depiction_lab_ui.py tests/test_depiction_lab_export.py`: **41 passed** (1.95 s; rerun final tras los cambios del reporte/perfiles: 41 passed en 1.95 s).
- Incluye 15 parámetros que cambian la salida/propiedad de un `BondItem` real; invariantes de grafos; independencia Current/Working; 12 escenas deterministas; overlays descartables; esquema JSON, migración v0→v1, carga/guardado; UI sincronizada y apagado Qt en subproceso acotado; exportación PNG/SVG/JSON/Markdown y hoja Current/Working de los 12 casos.
- Exportación completa y ciclo UI se prueban con Qt offscreen en procesos aislados y timeouts externos. Ninguna aprobación estética o preset editorial se da por establecida.
- `python -m compileall src tests tools packaging`: **PASS**.
- Ruff focal con `F401,F811,F821,E722,E741`: **PASS**. Ruff global exacto conserva únicamente el F401 histórico `import math` sin uso en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; archivo fuera de alcance y sin cambios.
- `openspec validate add-molecular-depiction-visual-lab --strict`: **PASS**.
- `git diff --check`: **PASS**.

## Suite existente por shards

Se generó un plan de 8 shards con **2,326 node IDs únicos**, exactamente los del manifiesto revisado `tests/ci/expected_pytest_nodeids.txt` (41 IDs nuevos del laboratorio). La herramienta verificó cobertura JUnit exacta, todas las asignaciones, artefactos, timeouts y señales; no faltan tests. Resultados:

| Índice | Resultado | Duración |
|---:|---|---:|
| 0 | 291 passed | 46.9 s |
| 1 | 287 passed, 4 skipped | 53.6 s |
| 2 | 289 passed, 2 skipped | 46.8 s |
| 3 | 285 passed, 5 skipped, **1 fallo** | 274.7 s |
| 4 | 289 passed, 2 skipped | 138.4 s |
| 5 | 285 passed, 6 skipped | 267.5 s |
| 6 | 290 passed | 34.0 s |
| 7 | 290 passed | 17.0 s |

El único fallo de shard es `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`. En ejecución aislada en proceso fresco pasa (**1 passed, 0.22 s**). Es el mismo fallo Qt dependiente del orden que el baseline de esta campaña ya registraba en el shard 4; con el manifiesto ampliado la repartición lo ubica en shard 3. No lo causó el laboratorio, no se editó ni se ocultó. El plan completa 2,306 passed + 19 skipped + 1 fallo = 2,326; el verificador reporta esa excepción sin error de cobertura. Los ocho shards sí se ejecutaron; no se repitió `pytest -q` monolítico después del baseline: el run inicial se agotó al límite externo de 600 s cerca de 34 %.

Un intento de herramienta con índice 8 fue rechazado por la validación del rango (shards indexados 0–7); no inició tests. Se corrió índice 0 a continuación y los ocho resultados quedaron completos.

## Límites y revisión

- El resultado del fallo de shard deja la suite global con un fallo heredado; la aceptación focal del laboratorio pasa. No se afirma CI verde completo.
- El renderer no usa un servicio orbital. En la base sólo se encontró `tools/visual/verify_orbitals_rendering.py`, una utilidad de verificación, no una GUI de calibración reusable; por ello no se agregó calibración de orbitales.
- B4 se identifica como `geometric-fixture`, no como molécula. No se ejecuta Clean2D.
- No existe aprobación visual, elección de estilo editorial ni promoción a defaults de producción. Evidencia estática y reporte de procedencia sirven para revisión posterior del propietario.
