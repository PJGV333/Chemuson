# Tareas — SIGABRT Qt/QUndoStack

## A. Baseline y reproducción
- [x] A.1 Registrar SHA actual, estado limpio, comandos baseline y el fallo shard 6 de Actions/local.
- [x] A.2 Construir snapshots base `61f7647b` y actual con `git archive`, sin worktrees; verificar import paths y mismo runtime.
- [x] A.3 Ejecutar la misma lista/orden de IDs del shard 6 en base y actual, omitiendo sólo el test nuevo que no existe en base; guardar logs y exits.
- [x] A.4 Ejecutar test aislado, archivo UI aislado y comparar el prefijo/orden que precede al abort.

## B. Causa y fix
- [x] B.1 Auditar el candidato `CanvasTabManager.cleanChanged` y el lifecycle observado; el trigger reproducible se ubica en el flush global de DeferredDelete del test.
- [x] B.2 Clasificar evidencia como defecto de aislamiento del test + orden/lifecycle Qt, preexistente en ambos SHAs; no regresión de producción ni sólo entorno.
- [x] B.3 Cambiar sólo el receiver del `DeferredDelete` al overlay; sin cambios de producción, skips/xfails/timeouts ni planner.
- [x] B.4 Confirmar que el cambio sólo modifica el cuerpo de un test existente: colección sigue en 2.598 IDs, el manifiesto no cambia y `tests/test_test_workflow_contract.py` pasa.

## C. Cierre
- [x] C.1 Focal UI/workflow 42 passed; A/B ordenado green en ambos snapshots; shard 6 actual completo PASS.
- [x] C.2 Ejecutar los ocho shards secuencialmente; verifier confirma JUnit/timings/IDs exactos.
- [ ] C.3 Verificar Windows smoke y Flatpak smoke en Actions sobre SHA final.
- [ ] C.4 Registrar evidencia y dictamen; suite completa verde es condición de cierre.
- [x] C.5 No hacer merge, tag, Release ni publicación; trabajar sólo en `stabilize/v0.3.0-beta.1`.
