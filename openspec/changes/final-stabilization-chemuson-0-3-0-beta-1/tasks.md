# Tareas — estabilización final ChemUSON 0.3.0-beta.1

## A. Auditoría y estabilización CI
- [x] A.1 Confirmar rama base/SHA/estado, leer AGENTS, arquitectura, OpenSpecs y documentación de publicación; crear `stabilize/v0.3.0-beta.1` desde la beta esperada en el mismo directorio.
- [x] A.2 Capturar baseline antes de editar y registrar limitaciones del intérprete local.
- [x] A.3 Descargar plan, JUnit, logs, timings y reportes de los ocho shards del run `39679607676`; reconciliar conteos y fallos exactos.
- [x] A.4 Actualizar el doble `HeldCompChemWorker` con `progress_ready`; el test de cierre alcanza y pasa sus aserciones.
- [x] A.5 Diagnosticar y corregir el acoplamiento de `thread.quit` al event loop GUI con una regresión red/green; conservar el timeout de 3 s. Confirmado en Actions del SHA publicado `5ca26395`.
- [x] A.6 Ejecutar pruebas focalizadas y los ocho shards. Actions #39782991976 (CI #422 informado por el owner) verificó el SHA exacto `f0151cf09a7200acc53cb6b371105e79fd3cac1d`: 2.598/2.598 IDs, 2.576 passed, 22 skipped, cero fallos/errores y los ocho JUnit/timings completos. Riesgo persistente: shard 1 duró 299.398 s del límite de 300; el run anterior #39724355063 sobre código funcionalmente idéntico había agotado el límite en shards 5/7. Esto no acredita estabilidad estadística. Detalle en `validation.md`.
- [x] A.7 OpenSpec estricto, arquitectura incluida en la matriz (283 tests), compileall, Ruff (sólo F401 basal) y `git diff --check` verificados.
- [ ] A.8 Recalibrar el perfil con evidencia de Actions reciente. El perfil `d2d7db19…5341e` cubre 253 archivos/fallback 0.2 s y pasó una matriz remota (#39723284720), pero #39724355063 mostró dos timeouts con la misma asignación. El candidato experimental `fdb58a09…3254e`, basado en timings completos/parciales de esos runs, predice 268.74 s equilibrados; su matriz local tuvo timeout en shard 2 y no se adopta. Falta diagnóstico/validación sin cambiar el cap de 300 s ni alterar el manifiesto.

## B. Auditoría de distribución
- [x] B.1 Revisar workflows/scripts/validadores existentes y documentar la diferencia entre el Build Preview Type 2 auténtico y las menciones históricas al ejecutable renombrado.
- [ ] B.2 Validar un Build Preview nuevo de Windows portable, Inno Setup, AppImage Type 2 y Flatpak para el mismo SHA candidato. El run #39679607716 generó los cuatro formatos, pero su SHA `61f7647b…` no es el actual `f0151cf…`; no se acepta como set del candidato.
- [ ] B.3 La auditoría por formato de #39679607716 verificó hashes/provenance, Type 2/extracción, 69 SVG, ChemName/RDKit/3D y recursos Flatpak. Se mantiene pendiente para el candidato exacto porque entre ambos SHAs cambió `src/chemuson/gui/controllers/compchem3d_controller.py`; Inno no se instaló y Flatpak sólo se importó/inspeccionó en repositorio OSTree temporal.
- [x] B.4 Registrar el set exacto `f0151cf…` como **NOT BUILT** y las verificaciones locales/visuales inaccesibles como **NOT TESTED**, sin inferir PASS de otro SHA.

## C. Aceptación química y gráfica
- [x] C.1 Mantener sin cambios la matriz existente `docs/release/manual-acceptance-0.3.0.md`; no elevar resultados manuales ni prioridades sin evidencia identificable.
- [x] C.2 La matriz automatizada incluida en CI #39782991976 pasó: 417/417 tests relacionados con Clean2D; 19/19 amida/urea, 14/14 MOL P0, 12/12 importación estereoquímica SMILES, 13/13 Live MM stream, 5/5 proyección estereoquímica, 58/59 CompChem 3D, 35/36 geometry service y 71/76 casos ChemName acceptance. Los skips y límites quedan desglosados en `validation.md`; no equivale a aceptación manual universal.
- [x] C.3 Separar smoke automatizado, inspección estática y aceptación manual. La prueba Qt/OpenGL real, Windows interactivo, instalación y validación gráfica de los paquetes exactos se registran NOT TESTED.
- [ ] C.4 Ejecutar aceptación de los cuatro paquetes exactos sólo cuando existan artifacts con el SHA candidato y el owner pueda probarlos; mantener P0/P1 pendientes.

## D. Preparación de publicación
- [ ] D.1 Se actualizan notas, OpenSpec y guía preview con SHAs, runs, artifact IDs y hashes. La matriz `manual-acceptance-0.3.0.md` queda deliberadamente intacta por instrucción del owner; por ello esta tarea documental completa permanece abierta.
- [x] D.2 Dictamen final: **NO-GO** para publicación; CI verde registrado, pero faltan builds exactos `f0151cf…`, aceptación manual P0/P1 y un margen de CI confiable. Advertencias y bloqueos en `validation.md`.
- [x] D.3 No se hizo merge, worktree, tag, Release, PR ni cambio a canales públicos/`main`/`gh-pages`; cualquier actualización se limita a la rama autorizada de estabilización.
