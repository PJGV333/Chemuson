# Validación focalizada — editor molecular 3D interactivo

Fecha: 2026-10-16. Rama: `feature/chemuson-3d`. HEAD previo a estos commits: `56328663dfc7b54ffaba7f9d7d60bab2dbf067c7`. Se preservaron los 18 paths que ya estaban modificados/no rastreados; este informe se añade como documentación de evidencia. No se hizo push, merge ni PR.

## Investigación de fase 3

- Se mantiene `avogadro-1.2.0-analysis.md`, limitado a inspección estática identificada por tag/commit, lectura del código fuente y observaciones de Ghidra. Distingue hechos de fuente, inferencias estáticas y pseudocódigo; no afirma ejecución de Avogadro ni incluye código GPL.
- La evidencia de baseline continúa en `baseline.md` y `baseline/`. El `pytest -q` baseline terminó con SIGSEGV alrededor del 61%; esta sesión no repitió la suite completa.
- El diseño mantiene `MolGraph` como autoridad química, coordenadas 3D por canvas/IDs estables, H generados efímeros y previews fuera de persistencia/undo.

## Verificaciones automatizadas ejecutadas

| Verificación | Resultado |
|---|---|
| `pytest -q tests/test_geometry3d_editing_metrics.py tests/test_geometry3d_commands.py tests/test_geometry3d_persistence.py tests/test_geometry3d_service.py tests/test_viewport3d_widget.py tests/test_viewport3d_camera.py tests/test_viewport3d_picking.py tests/test_compchem3d_dock.py tests/test_rdkit_live_stream.py tests/test_gui_async_worker_shutdown.py` | 74 passed, 1 skipped, 12.98 s |
| `pytest -q tests/test_geometry3d_persistence.py tests/test_canvas_opacity.py tests/test_main_window_tabs.py` | 22 passed, 6.74 s |
| `pytest -q tests/architecture/test_import_boundaries.py tests/architecture/test_module_catalog.py` | 84 passed, 6.37 s |
| `ruff check` sobre los archivos Python del cambio con `--select F401,F811,F821,E722,E741` | All checks passed |
| `python -m compileall -q` sobre archivos Python afectados | exit 0 |
| `openspec validate add-editable-3d-molecular-editor --strict` | Change valid |
| `git diff --check` | exit 0 |

No se suman los conteos de pytest: los grupos se solapan. El Ruff selectivo sigue el gate del repositorio. Un intento distinto de Ruff sin `--select` incluyó `architecture/modules.yml` como si fuera Python y reglas ajenas al gate; produjo miles de diagnósticos no accionables y no se usa como verificación del cambio.

## Demostración local de butano

En la primera ejecución gráfica en `DISPLAY=:0`, la ventana ChemUSON real inicializó el contexto OpenGL (`gl_ready: true`) y mostró butano con 14 átomos/13 enlaces (C4H10, con H derivados de RDKit). Se usaron los controles visibles del dock para aplicar al cuarteto C–C–C–C una torsión de 60° y habilitar `Live MM (RDKit)` con MMFF94. La activación se produjo desde los widgets Qt del dock; el arrastre con ratón se verifica separadamente en pruebas de widget.

La secuencia real del worker aislado NDJSON entregó 30 frames distintos, de los pasos 8 a 240. En todos los frames el controlador aún tenía activo el job al entregar el evento a UI. Las coordenadas cambiaron en los 30 snapshots; la energía descendió en 29 transiciones consecutivas, de 281.439545887 a −4.293762924 kcal/mol. El dihedro observado fue 133.239294° inicial, 60.000000° tras la operación del dock y 65.287969° al terminar la relajación. El resultado final se confirmó en el canvas y quedó como una operación undoable.

Evidencia no versionada (scratch local):

- Captura visual limpia del viewport y los controles: `/home/ccachyavgp/.hermes/cache/scratch/chemuson-phase4-butane.png`.
- Secuencia completa de energías/coordenadas: `/home/ccachyavgp/.hermes/cache/scratch/chemuson-phase4-butane.json`.

### Límite de cierre del proceso

La ejecución que produjo esa captura y secuencia imprimió el resultado funcional, pero terminó con exit 139 (SIGSEGV) durante el desmontaje de wrappers Qt/SIP. `coredumpctl info` para PID 1490286 mostró `sip_api_get_address` → `cleanup_qobject` → `sip_api_visit_wrappers` → `release_QApplication`; no se estableció causa raíz ni se atribuye a esta implementación. Es un traceback distinto del SIGSEGV del baseline pytest (que ocurre en `threading._worker`). Una repetición con el event loop habitual fue bloqueada por la aprobación adicional de la herramienta; no se intentó otra ruta ni se volverá a lanzar desde el agente.

Por esa razón, la interacción y la relajación en vivo están demostradas, pero no se certifica un cierre de proceso limpio ni se marca completa la validación gráfica. Etano y benceno tampoco se validaron gráficamente en esta sesión. La demostración manual debe confirmar esas cuestiones.

## Pasos mínimos de verificación manual en CachyOS

Desde la raíz del repositorio, ejecuta la prueba local ya preparada (sin privilegios y sin escribir configuración fuera de scratch):

```bash
env DISPLAY=:0 QT_QPA_PLATFORM=xcb \
  XDG_CONFIG_HOME="$HOME/.hermes/cache/scratch/chemuson-phase4-xdg-config" \
  XDG_CACHE_HOME="$HOME/.hermes/cache/scratch/chemuson-phase4-xdg-cache" \
  XDG_DATA_HOME="$HOME/.hermes/cache/scratch/chemuson-phase4-xdg-data" \
  CHEMUSON_DEMO_SCREENSHOT="$HOME/.hermes/cache/scratch/chemuson-phase4-butane.png" \
  PYTHONPATH=src ./.venv/bin/python \
  "$HOME/.hermes/cache/scratch/chemuson_phase4_live_gui_demo.py"
```

En la UI, comprueba que el dihedro llegue a 60°, activa MMFF94 Live MM, observa varios pasos con energía/coordenadas cambiantes antes de `done`, y verifica undo/redo de torsión y geometría final. Revisa que la aplicación cierre sin SIGSEGV. Para completar la validación gráfica, abre además etano y benceno y confirma que renderizan sin modificar `MolGraph` con hidrógenos derivados.

## Pendiente explícito

- Repetir manualmente el cierre gráfico; el intento posterior del agente quedó bloqueado y no debe reintentarse mediante otra herramienta.
- Smoke gráfico de etano y benceno.
- No se ejecutó la suite completa por el SIGSEGV de baseline conocido.

## Consolidación y pulido de fases 1–4 — addendum 2026-10-16

Este addendum amplía la evidencia sin editar ni sustituir los resultados de investigación, baseline o demostración gráfica anteriores. El trabajo se realizó en `feature/chemuson-3d`, HEAD inicial de esta continuación `a7dc0dcadf3d6a5fe735ed556f23f1ab5373cf6e`, upstream `origin/feature/chemuson-3d`. No hubo push, merge ni PR.

### 1. Defectos corregidos

- El constructor 3D valida la valencia sobre una copia candidata antes de abrir la macro de undo. El rechazo comunica el motivo y no altera el modelo, las coordenadas ni el historial.
- `AddBondCommand` rechaza un átomo inicial/final inexistente, autoenlace o par duplicado antes de modificar el grafo. La prueba existente confirma que agregar átomo+enlace sigue siendo una sola acción reversible.
- Las regresiones de valencia incompatible, endpoint inexistente y duplicado comparan el estado de grafo/geometría/historial antes y después.

### 2. Interfaz y renderizado

No se cambiaron shaders, iluminación, paleta, tamaños, layout del dock ni navegación. La instrucción de validación gráfica anterior prohíbe reintentar desde el agente la demostración bloqueada; por ello no presento cambios visuales ni una nueva verificación OpenGL como hechos. La mejora visual y la revisión de compresión/redimensionamiento del dock permanecen pendientes de una prueba gráfica autorizada.

### 3. Estado de integración 2D ↔ 3D

- El control existente `Generar conformero 3D` y el worker RDKit aislado siguen usando el `MolGraph` documental; se probó el servicio sobre grafos 2D importados desde SMILES. El generador produce coordenadas de conformero RDKit y una escena con H derivados sin añadirlos al grafo.
- La construcción directa conserva `MolGraph` como autoridad. La prueba crea átomo+enlace en un único paso undoable; las coordenadas 3D no mueven el layout 2D. La proyección 3D→2D requiere confirmación explícita y se puede deshacer.
- La geometría `.cmsn` se guardó y reabrió; IDs documentales y layout 2D se conservaron, mientras los H derivados no se serializaron.
- Se probó una estructura con carga formal +1 por el servicio 3D: `MolGraph` mantuvo la carga y los H generados siguieron fuera del documento.
- **Límite estereoquímico:** estos tests no validan que la conformación 3D respete cuñas/CIP. La inspección del payload confirma que `BondStereo` wedge/hash no se serializa al worker y `stereo_cip` no se convierte en una etiqueta quiral RDKit; por tanto no se afirma preservación de la geometría estereoquímica. El grafo documental no se muta, pero esta condición sigue abierta.

### 4. Matriz química ejecutada

Los cinco casos se importaron como `MolGraph`, se generó conformero y se comparó la escena con el documento. Todos los valores de escena fueron finitos; los IDs no generados coincidieron con los IDs del `MolGraph`.

| Caso | Átomos/enlaces documentales | Átomos/enlaces de escena | Resultado adicional |
|---|---:|---:|---|
| Agua | 1 / 0 | 3 / 2 | 2 H derivados; ángulo H–O–H observado 103.98° |
| Etano | 2 / 1 | 8 / 7 | 6 H derivados |
| Etanol | 3 / 2 | 9 / 8 | 6 H derivados; C–C/C–O documentales y O–H derivado |
| Butano | 4 / 3 | 14 / 13 | 10 H derivados; torsión/Live MM cubierta por pruebas focalizadas |
| Benceno | 6 / 6 | 12 / 12 | 6 H derivados; anillo cíclico 3 enlaces sencillos + 3 dobles (Kekulé) |

Rangos Z medidos en la escena (Å): agua 0.348, etano 2.579, etanol 2.316, butano 3.395 y benceno 4.710. Ninguna escena fue una proyección con todos los Z iguales a cero.

En la ruta de importación probada, benceno está representado como Kekulé alternante; `is_aromatic` no está marcado en sus seis enlaces. No se interpreta esa salida como preservación de flags aromáticos explícitos.

Live MM ejecutó MMFF94, MMFF94s y UFF sobre butano: cada flujo emitió 40 snapshots (pasos 2–80), con coordenadas/energía progresivas y átomo fijado. Energías del primer al último frame en esa corrida: MMFF94 3917.374218→−2.605608; MMFF94s 3917.374218→−2.605608; UFF 1154.252171→5.212458 kcal/mol. El test del worker verifica snapshots distintos, finitud, paso/energía, restricción del átomo y resultado terminal.

### Seguimiento: continuidad de hidrógenos derivados

La revisión de los callbacks de relajación detectó que las escenas de Live MM reconstruían cada frame sin los metadatos de H derivados. Se corrigió el worker/backend para incluir y conservar posiciones y enlaces de esos H en cada frame NDJSON; las vistas previas de movimiento y las restricciones geométricas ahora trasladan cada H con su átomo documental padre. Pruebas headless verifican que todas las escenas de butano en MMFF94/MMFF94s/UFF conserven 14 átomos y 13 enlaces, y que el preview de movimiento mantenga H unido. Los H siguen siendo metadata efímera, fuera de `MolGraph` y del archivo persistido.

Tras cambios de topología, los H previos no se reutilizan porque su valencia/adjacencia puede quedar obsoleta; reaparecen al generar o relajar la geometría con RDKit. No se verificó visualmente OpenGL en esta actualización.

### 5. Rendimiento observado

Una corrida de conformer RDKit con proceso aislado midió 0.164–0.239 s por molécula para estos casos pequeños; construir `SceneMolecule3D` midió 0.000067–0.000106 s. Esto apunta a generación/proceso RDKit, no a ensamblaje de escena, como el costo dominante de esta muestra. No se midieron FPS, interacción OpenGL, costo adicional de serializar H en frames ni moléculas grandes; no se aplicaron optimizaciones de renderizado sin esa evidencia.

### 6. SIGSEGV y cierre

El shard `tests/test_gui_async_worker_shutdown.py` pasó **2 tests** en modo Qt offscreen y terminó con exit 0; el shard combinado focalizado también terminó limpiamente. Esto no reproduce el cierre de la demostración real con contexto OpenGL. La causa registrada anteriormente (`sip_api_get_address`/limpieza de wrappers) permanece sin determinar y distinta del fallo de worker `threading._worker` de la suite baseline. No se ejecutó otra demo gráfica ni se repitió pytest completo: el informe/tareas anteriores documentan la aprobación bloqueada y el límite de suite, respectivamente. El SIGSEGV no se declara resuelto.

### 7. Validaciones, limitaciones y pasos para CachyOS

- Suite focalizada de viewport, geometría, persistencia, comandos, Live MM, pestañas y shutdown: **102 passed in 43.76s**; no equivale a la suite completa.
- Tests arquitectónicos: **280 passed in 21.62s**.
- `compileall`: exit 0; Ruff scoped a los archivos Python modificados: All checks passed; `git diff --check`: exit 0.
- Ambos `openspec validate ... --strict` indicaron `Change valid`. El CLI emitió avisos de formato `rules` en artefactos y los ignoró; no se modificaron retrospectivamente.
- Siguen abiertos: demostración gráfica completa y cierre limpio Qt/SIP; pulido visual/viewport; suite monolítica; preservación de geometría quiral; flags aromáticos explícitos en la ruta de importación probada; recorridos manuales completos de guardado/pestañas.
- Prueba manual mínima en CachyOS, desde la raíz del repositorio y con dependencias instaladas: `QT_QPA_PLATFORM=xcb PYTHONPATH=src ./.venv/bin/python -m chemuson`. Después, verifica los dos recorridos, guardar/reabrir, cambio de pestañas, Live MM y cierre limpio. Estas comprobaciones manuales no se consideran realizadas por esta ejecución.

Commits locales de código/tests: `58383965` (`fix(3d): reject invalid molecular builder edits`), `1b7fb1ec` (`test(3d): add chemistry and live-mm acceptance coverage`) y `4014d4f7` (`fix(3d): keep derived hydrogens in live frames`). Esta actualización de OpenSpec se registra en un commit documental separado. No se publica la rama.
