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

- Pruebas negativas dedicadas del constructor 3D: endpoint inexistente, enlace duplicado y rechazo de valencia incompatible.
- Repetir manualmente el cierre gráfico; el intento posterior del agente quedó bloqueado y no debe reintentarse mediante otra herramienta.
- Smoke gráfico de etano y benceno.
- No se ejecutó la suite completa por el SIGSEGV de baseline conocido.
