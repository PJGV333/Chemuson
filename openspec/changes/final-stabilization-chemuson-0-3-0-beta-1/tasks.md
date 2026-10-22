# Tareas — estabilización final ChemUSON 0.3.0-beta.1

## A. Auditoría y estabilización CI
- [x] A.1 Confirmar rama base/SHA/estado, leer AGENTS, arquitectura, OpenSpecs y documentación de publicación; crear `stabilize/v0.3.0-beta.1` desde la beta esperada en el mismo directorio.
- [x] A.2 Capturar baseline antes de editar y registrar limitaciones del intérprete local.
- [x] A.3 Descargar plan, JUnit, logs, timings y reportes de los ocho shards del run `39679607676`; reconciliar conteos y fallos exactos.
- [x] A.4 Actualizar el doble `HeldCompChemWorker` con `progress_ready`; el test de cierre alcanza y pasa sus aserciones.
- [x] A.5 Diagnosticar y corregir el acoplamiento de `thread.quit` al event loop GUI con una regresión red/green; conservar el timeout de 3 s. Confirmado en Actions del SHA publicado `5ca26395`.
- [x] A.6 Ejecutar pruebas focalizadas y los ocho shards: el SIGABRT anterior de shard 6 y el recién expuesto de shard 5 se reprodujeron y corrigieron como aislamiento de DeferredDelete en tests; matriz local final 2.598 IDs, 2.576 passed, 22 skipped, cero fallos, JUnit/timings exactos, máximo 177.566 s. Actions #39723284720 sobre `5ca26395` pasó plan, los ocho shards, resumen y smoke Windows/Flatpak; 2.598/2.598 JUnit, máximo 293.063 s (<300 s).
- [x] A.7 OpenSpec estricto, arquitectura incluida en la matriz (283 tests), compileall, Ruff (sólo F401 basal) y `git diff --check` verificados.
- [x] A.8 Actualizar el perfil de duraciones con mediciones trazables de Actions #39697199469/#39679607676/#38996089270; comprobar 253 archivos, fallback 0.2 s/test y plan exacto determinista. Ocho shards locales confirman todos los tiempos <300 s; no se amplió timeout ni se cambió manifiesto.

## B. Auditoría de distribución
- [ ] B.1 Revisar workflows/scripts/validadores existentes y contradicciones en documentación histórica AppImage.
- [ ] B.2 Validar evidencia Build Preview nueva para Windows portable, Inno Setup, AppImage Type 2 y Flatpak, todos con mismo SHA; no publicar.
- [ ] B.3 Comprobar versión/identidad/recursos SVG y ChemName/RDKit/3D, AppRun, extracción sin FUSE, aislamiento Flatpak, checksum y provenance; registrar límites de host.
- [ ] B.4 Registrar builds inexistentes o verificaciones inaccesibles como NOT BUILT/NOT TESTED, no PASS inferido.

## C. Aceptación química y gráfica
- [ ] C.1 Mantener la matriz existente `docs/release/manual-acceptance-0.3.0.md` y priorizar inicio/cierre/workers, persistencia/autosave/recovery y casos de identidad/estereoquímica.
- [ ] C.2 Verificar dibujo 2D, Clean2D/amidas/ureas/aromáticos, undo/redo, proyección 2D→3D/tetrandrina y Live MM con evidencia disponible.
- [ ] C.3 Separar los smoke automatizados de aceptación manual; Windows/OpenGL real que no pueda probarse queda NOT TESTED.
- [ ] C.4 Ejecutar aceptación de cuatro paquetes reales sólo cuando el owner tenga artifacts identificables; conservar P0/P1 pendientes.

## D. Preparación de publicación
- [ ] D.1 Actualizar notas, matriz, OpenSpec y guías de instalación con SHA/run/artifact/checksums exactos.
- [ ] D.2 Presentar dictamen GO/NO-GO y lista P0/P1, advertencias aceptables, resultados CI y pruebas pendientes.
- [ ] D.3 No hacer merge, crear worktrees, tag, Release, PR ni cambiar canales públicos; commit y push normal sólo a la rama de estabilización con autorización explícita del owner (autorizado en esta continuación).
