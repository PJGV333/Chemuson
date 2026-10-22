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
- [x] C.3 En Actions sobre SHA final, `windows-smoke` y `flatpak-smoke` PASS; el workflow global sigue fallido por timeouts en shards 5 y 6.
- [x] C.4 Registrar evidencia y dictamen **NO-GO**: el SIGABRT queda corregido, pero la matriz remota completa no alcanzó cobertura/JUnit por los timeouts.
- [x] C.5 No hacer merge, tag, Release ni publicación; trabajar sólo en `stabilize/v0.3.0-beta.1`.

## D. Incidente adicional de shard 5 — ciclo de vida Molecular Assistant
- [x] D.1 Capturar el estado Git inicial y respaldar el diff fuera del repositorio antes de editar; conservar `.tmp-beta1-readiness/` durante la investigación.
- [x] D.2 Reproducir el shard 5 de 349 IDs en el mismo orden sobre snapshot base `3f800c7`: SIGABRT tras 272 timings, sin JUnit.
- [x] D.3 Aplicar en snapshot la entrega inicial del evento al QObject objetivo; el mismo shard termina 349/349, 341 passed, 8 skipped, JUnit/timings completos.
- [x] D.4 Aplicar receiver-scoping a helpers de lifecycle y probes Qt equivalentes, sin tocar producción y manteniendo las aserciones.
- [x] D.5 Ejecutar focales y matriz local completa: 2.598 IDs, 2.576 passed, 22 skipped, cero fallos; ocho JUnit/timings completos, máximo 177.566 s (<300 s).
- [ ] D.6 Ejecutar y registrar Actions sobre el SHA publicado; los resultados remotos previos siguen fallidos/incompletos hasta un nuevo run.
