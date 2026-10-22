# Diseño — diagnóstico causal y comparación aislada

## Causa clasificada

El trigger es del **harness + orden/lifecycle Qt**, no una regresión de producción en `CanvasTabManager` ni una incompatibilidad aislada de runner. `test_onboarding_close_releases_overlay` llama `QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)`, que vacía los DeferredDelete de toda la aplicación mientras sólo intenta destruir su overlay. En el shard ordenado, eso entrega deletes pendientes de otros widgets/undo stacks y Qt aborta dentro de `QUndoStack::~QUndoStack`/`cleanChanged`. El archivo UI aislado no tiene el mismo backlog.

La comparación en snapshots limpios de base y candidato, mismo Python 3.11.15/PyQt6 6.11.0/Qt 6.11.2/RDKit 2026.9.1, mismo cwd aislado y los mismos 372 IDs, da SIGABRT en la misma línea en ambos. Un cambio experimental en **sólo** ese test (`sendPostedEvents(overlay, DeferredDelete)`) hace que ambos snapshots terminen 372/372, 370 passed, 2 skipped y JUnit completo. No hace falta tocar la lambda de producción; el cambio evita que el test descargue eventos de otros owners.

## Hipótesis ya descartadas

1. No se demostró una regresión de `aff565a` frente a `61f7647b` para el abort de onboarding: el test aborta con la misma secuencia en ambos.
2. No es sólo un problema de Python/Qt local: el Actions run `39689871971` reproduce el mismo test con Python 3.11.17; el A/B local con Python 3.11.15/Qt 6.11.2 conserva el mismo resultado base/candidato.
3. El `DirectConnection` CompChem no es la causa del abort de onboarding: los dos árboles, con la misma secuencia, fallan igual antes del cambio de receiver. La ruta CompChem pasa en ambos durante el A/B Python 3.11.
4. El cambio de receiver se comprobó experimentalmente en snapshots temporales; aún debe entrar al checkout, ejecutarse focal/shards y verificarse en Actions.

## Experimento A/B

- Crear snapshots fuente con `git archive` de `61f7647b` y `aff565a064f17dccf98b0b4f23c7ddeba593ca31`, fuera del índice y sin `git worktree`; confirmar que `chemuson.__file__` de cada proceso apunta al snapshot correspondiente.
- Usar el mismo intérprete/dependencias/QPA y el mismo orden de IDs del shard 6 actual. Excluir únicamente del set compartido el nuevo node ID de regresión, que no existe en el base; no excluir ningún test preexistente.
- Ejecutar cada snapshot bajo timeout externo 300 s y preservar exit status, log y el último test. Repetir sólo donde aporte una comparación causal; el timeout nunca convierte un fallo en PASS.
- Comparar además el test de onboarding aislado y el archivo `test_ui_polish.py` para separar fallo local del orden previo. Revisar conexiones, owners y eventos DeferredDelete antes de decidir un fix.

## Corrección mínima autorizada

Editar sólo `tests/test_ui_polish.py`: enviar el DeferredDelete al `overlay` exacto, no a `None`, y documentar que no se deben procesar deletes ajenos en esta prueba. Así se conserva la afirmación que prueba (el overlay se destruye) sin consumir el backlog de otros widgets. No cambiar `src/chemuson/gui/tab_manager.py`, el planner, asserts funcionales, skips ni timeouts.

## Diagnóstico añadido — shard 5, cierre de Molecular Assistant

La secuencia de shard 5 (349 IDs, perfil candidato) se ejecutó en snapshots aislados de `git archive` de HEAD `3f800c715ee937d9ac72385d8189bb9d96902ea3` y una variante con receiver específico. Se confirmó `chemuson.__file__` dentro de cada snapshot y se mantuvieron el mismo runtime Python 3.14.8 / pytest 9.1.1 / PyQt6 6.11.0 / Qt 6.11.2 / RDKit 2026.03.6 y `QT_QPA_PLATFORM=offscreen`.

- Base, helper original: exit -6/SIGABRT tras 272 timings; el test iniciado en índice 272 fue `test_unstarted_dialog_close_and_window_close_leave_no_registry_entries`; JUnit ausente, 93.597 s. El stack nativo atraviesa `_process_qt_deletes`, `QUndoStack::~QUndoStack` y `cleanChanged`.
- Receiver-scoped: el mismo shard/orden termina 341 passed, 8 skipped, 349/349 JUnit y timings, exit 0 en 151.017 s. En el checkout se repite el shard completo a 151.145 s; el archivo lifecycle pasa 9/9.
- Causa confirmada para este aborto: `_process_qt_deletes()` vaciaba con receiver `None` el DeferredDelete global al verificar una destrucción local, forzando la destrucción de QUndoStack/árboles QWidget no poseídos por el diálogo bajo prueba. La corrección entrega el evento directamente al diálogo/canvas/ventana objetivo y después procesa eventos ordinarios. No cambia producción ni aserciones de destrucción/registros/trabajadores.
- Un primer experimento pasó `QCoreApplication.processEvents()` antes de usar el receiver y produjo `RuntimeError: wrapped C/C++ object ... has been deleted`: el primer bombeo ya había destruido el diálogo. El helper final entrega primero el DeferredDelete al receiver y sólo después llama `processEvents()`.
- Las llamadas equivalentes en `tests/probes/p0_shutdown_workflow.py` y `tests/probes/p0_open_mol.py` pertenecen a procesos Qt secundarios independientes; se inspeccionaron y se limitaron también a widget/canvas/ventana exactos. Sus pruebas focales terminan 16/16. No se modificó código de aplicación ni otras suites indiscriminadamente.

## Verificación

Una corrección sólo se acepta si el test de regresión falla en el caso basal pertinente, pasa después, la secuencia previamente abortante termina con JUnit y timings completos, los ocho shards pasan con cobertura exacta, y Windows/Flatpak smoke pasan en Actions para el SHA final. Los resultados automáticos no son aceptación manual.
