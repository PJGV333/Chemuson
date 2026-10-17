# Evidencia de validación — estabilización 3D fases 1–4

## Registro del 2026-10-16

- Rama: `feature/chemuson-3d`; HEAD de partida: `a9937e7f245c9ce88ebd3f421fe2513db5b77203`; upstream: `origin/feature/chemuson-3d`; el branch local estaba 2 commits adelante. No se modificó `main` ni se publicó nada.
- Alcance: continuidad de H derivados, estereoquímica explícita, aromaticidad/órdenes, render/picking, estado de cierre y controles del viewport.
- No se lanzó la GUI ni se repitió la suite completa, por instrucción explícita. La aceptación visual y la causa del SIGSEGV siguen abiertas; esta evidencia no declara la fase 4 plenamente validada.

## Pruebas y verificaciones

- Regresiones focalizadas offscreen en geometría 3D, importación/estéreo, persistencia, Live MM, dock, viewport, apagado async y regresiones de canvas/tabs: **143 passed, 1 skipped, 1 deselected in 45.32s**. Se excluyó sólo `test_tetrandrine_import_preserves_visual_stereo` para obtener la ejecución verde.
- La misma selección con ese test incluido detectó su aserción fallida. El SMILES de Tetrandrine del fixture no contiene marcadores `@`; la ruta baseline de `rdkit_io.py` tomada de `HEAD` produjo **0 enlaces wedge/hash y 0 etiquetas estéreo**. Esto concuerda con no inventar configuración estereoquímica no especificada y muestra que el fallo no lo introduce el parser modificado. No se cambió el test ni se ocultó el fallo; requiere decisión/actualización del fixture fuera de esta corrección funcional.
- Pruebas de arquitectura (`test_import_boundaries.py`, `test_module_catalog.py`): **84 passed in 5.79s**.
- `python -m compileall -q src tests tools packaging`: **exit 0**.
- Ruff focalizado sobre módulos y pruebas modificados, reglas `F401,F811,F821,E722,E741`: **All checks passed**.
- OpenSpec `--strict` para `add-interactive-3d-molecular-viewport`, `add-editable-3d-molecular-editor` y este cambio: **los tres válidos**.
- `git diff --check`: **exit 0**, sin errores de whitespace.
- La suite completa sigue sin repetirse: existe el SIGSEGV histórico y el usuario pidió evitar otra ejecución innecesaria.

## Evidencia de continuidad y fidelidad química

Las regresiones focalizadas ejercitan transformaciones de geometría con H derivados ligados al átomo padre; cancelación de drag; controles de distancia/ángulo/dihedro; undo/redo; frames Live MM de butano; importación R/S y E/Z; rechazo de geometrías con configuración invertida; marcas estéreo ambiguas; hash químico diferenciado por estéreo; benzene aromático/Kekulé; persistencia `.cmsn`; carriles dobles/triples/aromáticos y picking sincronizado. El test de aceptación de moléculas incluye etano, etanol, butano y moléculas aromáticas/múltiples. Son pruebas automatizadas, no una verificación visual de la aplicación completa.

## Costo medido de malla (CPU)

Medición con `build_scene_mesh` después de generar conformeros RDKit, 1 calentamiento + 5 corridas por molécula. No abre Qt/OpenGL y no mide upload GPU ni FPS; por lo tanto no justifica cambios de buffers ni infraestructura gráfica.

| Molécula | Átomos | Enlaces | Vértices | Mediana CPU | Máximo de 5 corridas (p95 empírico) |
|---|---:|---:|---:|---:|---:|
| Etano | 8 | 7 | 10,104 | 4.176 ms | 4.769 ms |
| Butano | 14 | 13 | 17,736 | 7.465 ms | 8.381 ms |
| Benceno | 12 | 12 | 15,480 | 7.357 ms | 7.679 ms |
| Eicosano | 62 | 61 | 78,792 | 37.339 ms | 39.342 ms |

Script reproducible (solo construcción de malla): `/home/ccachyavgp/.hermes/cache/scratch/chemuson_mesh_cost.py`. La medición muestra un costo CPU mayor para la molécula de 62 átomos, pero no se optimizó la ruta: falta medir el ciclo de actualización/upload gráfico y validar visualmente bajo el flujo real de Live MM.

## Diagnóstico Qt/SIP: dos casos separados, causa no resuelta

1. **Cierre de pytest (core `476406`, 2026-10-15):** la pila nativa disponible sitúa el hilo principal en destrucción de `QGraphicsView`/`QObject` y `QUndoStack::~QUndoStack → clear → indexChanged → PyQtSlotProxy → PyQtSlot::call`. En ese snapshot, el worker Python aparece esperando en `queue.SimpleQueue.get`, no ejecutando RDKit. Esto descarta que ese snapshot demuestre un fallo directo dentro de una operación activa del worker; no descarta una relación indirecta ni identifica por qué se invoca el slot durante la destrucción.
2. **Demo gráfica 3D (core `1490286`, 2026-10-16):** termina con SIGSEGV durante liberación de wrappers SIP/QApplication; el frame superior registrado fue `sip_api_get_address` seguido por limpieza/visita de wrappers y liberación de QApplication. GDB no produjo backtrace adicional utilizable.

No hay evidencia suficiente para atribuir ninguno de los casos a buffers OpenGL, contexto gráfico, RDKit o una misma causa. Las pruebas offscreen de apagado async pasan, pero no reproducen ni descartan el fallo gráfico. No se relanzó la GUI para investigar.

### Procedimiento manual controlado para el equipo

No ejecutado por este agente. En un terminal gráfico interactivo, con el entorno del proyecto activo, registrar cada intento por separado:

1. Ejecutar y cerrar la aplicación sin abrir el viewport 3D; anotar si el proceso termina con código 0.
2. Ejecutar el caso reproducible de butano/Live MM que ya existe en el scratch, con `PYTHONFAULTHANDLER=1 DISPLAY=:0 PYTHONPATH=src ./.venv/bin/python /home/ccachyavgp/.hermes/cache/scratch/chemuson_phase4_live_gui_demo.py`; registrar código de salida y backtrace si falla.
3. Comparar cierre con viewport abierto pero sin trabajo, después de Live MM completado y —sólo en un intento separado— con trabajo activo. No confundir la reproducción del caso de demo con el SIGSEGV de pytest.
4. Para el caso pytest, usar `PYTHONFAULTHANDLER=1 QT_QPA_PLATFORM=offscreen ./.venv/bin/pytest -q` únicamente si el equipo decide repetir la suite conocida como inestable; capturar el nombre del último test y el core. El agente no volvió a ejecutarlo.

El script de demo provocó anteriormente un SIGSEGV al cerrar; correrlo puede terminar anormalmente. No se afirma que el fallo vaya a ser determinista. La causa raíz y la aceptación visual continúan pendientes.

## Cambios y commits locales

- `9cd14dc` — fidelidad estéreo de ChemIO/RDKit y clave química: `src/chemuson/chemio/_rdkit_worker.py`, `src/chemuson/chemio/rdkit_io.py`, `src/chemuson/chemio/rdkit_safe.py`, `src/chemuson/geometry3d/service.py`, `tests/test_geometry3d_service.py`, `tests/test_geometry3d_stereochemistry.py`, `tests/test_rdkit_io_stereo_import.py`.
- `13f2e03` — transforms de H derivados, semántica/enlaces de escena, mesh y picking: `src/chemuson/geometry3d/editing.py`, `src/chemuson/geometry3d/model.py`, `src/chemuson/gui/main_window.py`, `src/chemuson/gui/viewport3d/dock.py`, `src/chemuson/gui/viewport3d/mesh.py`, `src/chemuson/gui/viewport3d/picking.py`, `tests/test_compchem3d_dock.py`, `tests/test_geometry3d_aromaticity.py`, `tests/test_geometry3d_editing.py`, `tests/test_geometry3d_persistence.py`, `tests/test_rdkit_live_stream.py`, `tests/test_viewport3d_geometry.py`, `tests/test_viewport3d_picking.py`.
- `dbeda3d` — ciclo de vida del viewport: `src/chemuson/gui/viewport3d/widget.py`, `tests/test_viewport3d_widget.py`.
- El commit documental local incluye `AGENT_REPORT.md`, este `validation.md` y `tasks.md`. No se publicaron commits.