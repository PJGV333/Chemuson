# Validación — edición molecular 3D directa

**Fecha local:** 2026-10-16 (UTC−07:00) · **Rama:** `feature/chemuson-3d` · **HEAD de baseline:** `bfccb68c5e60221369f1ba83bfe701a4df3e8f8e`.

## Verificaciones ejecutadas

- `./.venv/bin/python -m compileall -q src tests tools packaging`: PASS.
- `./.venv/bin/pytest --collect-only -q`: PASS, **2201 tests collected in 0.47s**. Log: `/home/ccachyavgp/.hermes/cache/scratch/chemuson-campaign-final-collect.log`.
- Regresión GUI/editor: `tests/test_compchem3d_dock.py tests/test_viewport3d_widget.py`: **49 passed, 1 skipped** (la prueba RDKit de no duplicación se salta cuando el backend RDKit no está disponible en el entorno).
- Render/picking/cámara: `tests/test_viewport3d_geometry.py tests/test_viewport3d_picking.py tests/test_viewport3d_camera.py`: **28 passed**.
- Servicio, persistencia, edición geométrica y streaming: `tests/test_geometry3d_service.py tests/test_geometry3d_persistence.py tests/test_geometry3d_editing.py tests/test_rdkit_live_stream.py`: **44 passed, 1 skipped** (Open Babel no disponible; `tests/test_geometry3d_service.py:565`).
- Arquitectura: `tests/architecture`: **280 passed**.
- Ruff focalizado en todos los archivos de producción/prueba modificados, con `F401,F811,F821,E722,E741`: PASS (`All checks passed!`).
- `openspec validate improve-avogadro-inspired-3d-editing --strict`: PASS.
- `git diff --check`: PASS.

No se ejecutó `pytest -q` completo: el baseline previo de esta rama documenta SIGSEGV/exit 139 en la suite monolítica y el usuario pidió no repetirlo innecesariamente. El inventario de baseline y su salida están en `baseline.md`.

## Pruebas manuales y límites

- No se inició la GUI en esta campaña. La aceptación visual de estilos/iluminación y el recorrido en CachyOS siguen **pendientes**; ejecutar únicamente según `../stabilize-3d-phase-1-4-integrity/manual-acceptance-cachyos.md`, sin privilegios ni cambios de sistema. Esa guía también advierte el SIGSEGV histórico al cerrar; guardar primero el `.cmsn` temporal.
- No se ejecutaron binarios de Avogadro ni análisis dinámico/REA. La comparación se basa en inspección estática del árbol fuente local; la matriz detalla procedencia y límites.
- La verificación de no duplicación de H con RDKit no se ejecutó porque `_rdkit_available()` es falso en este entorno. La promoción, undo/redo, Save/Load y escena sin H duplicado se probaron por las rutas ChemUSON disponibles.
- No se midieron presupuestos de latencia monotónica del preview/streaming.

## Alcance funcional probado

Las regresiones cubren sincronización/debounce y descarte de geometría obsoleta; estado 3D por pestaña; edición de átomo/fragmento, construcción y cambios de orden de enlace undoables con validación de valencia; preview/drag/cancelación y promoción del H derivado; uso de la coordenada manipulada y restricción del H promovido en Live MM; medición directa de distancias/ángulos; estilos ball-and-stick, wireframe y space-filling con picking coherente; y el contrato estructural de iluminación en espacio de vista.

## Decisiones de arquitectura

- `MolGraph` conserva autoridad topológica/química; los H derivados usan ID visual no documental y solo se promueven con `AddAtomCommand` + `AddBondCommand` + `SetGeometry3DCommand` dentro de una macro undoable.
- La interacción 3D queda en el viewport/dock existentes; no se añadieron paquetes ni dependencias externas. La sincronización sigue filtrando resultados por snapshot químico y documento.
- Los estilos cambian malla y picking de forma conjunta; la iluminación se adapta a la vista. Los tests verifican contratos de datos/shader, no calidad visual.
- La comparación Avogadro es estática y no copia código GPL.

## Archivos modificados

- Producción: `src/chemuson/geometry3d/model.py`, `src/chemuson/gui/main_window.py`, `src/chemuson/gui/shell/assembly.py`, `src/chemuson/gui/viewport3d/{dock.py,mesh.py,picking.py,renderer.py,widget.py}`.
- Pruebas: `tests/test_compchem3d_dock.py`, `tests/test_viewport3d_geometry.py`, `tests/test_viewport3d_picking.py`, `tests/test_viewport3d_widget.py`.
- OpenSpec: `proposal.md`, `design.md`, `tasks.md`, `baseline.md`, `avogadro-source-matrix.md`, `roadmap-advanced-surfaces.md`, `validation.md` y las cinco especificaciones bajo `specs/`.

## Pendientes explícitos

Quedan abiertos los casos exactos butano→vacío→etanol y ráfagas vacías documentales; el paso RDKit de no duplicación; cambio del elemento de un átomo existente; handles directos para longitud/ángulo/torsión y selección visual guiada de restricciones; una prueba de orientación de cámara ligada a iluminación; y mediciones de latencia. No se presenta ninguna propuesta, documento de roadmap ni prueba textual de shader como aceptación gráfica o científica.
