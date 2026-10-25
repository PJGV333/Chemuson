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
- [x] B.2 Validar Build Preview nuevo para los cuatro paquetes, todos del SHA exacto `a38fa8045c858f8ff4b6002950f57e50292aea81`, rama remota `release/v0.3.0-beta.1-rc-prep`; Actions #39789295451 completó prepare, Windows portable/installer, AppImage Type 2, Flatpak y reporte con `success`, `publication=false`.
- [x] B.3 Verificar manifests, provenance/checksum de los cuatro paquetes, plataforma/versión, AppImage Type 2/extracción, 69 SVG, ChemName/RDKit/3D y Flatpak. Los smoke automatizados y las inspecciones están en `validation.md`; Inno no se instaló, ni se instaló Flatpak en el host. Esas pruebas de instalación/GUI permanecen NOT TESTED.
- [x] B.4 Registrar artifact IDs, SHA-256, provenance y resultados. Las cuatro descargas de Actions son los paquetes candidatos exactos; los artefactos históricos de `61f7647b…` no se usan como candidatos. Los entornos/instalaciones no ejecutados siguen **NOT TESTED**.

## C. Aceptación química y gráfica
- [x] C.1 Mantener sin cambios la matriz existente `docs/release/manual-acceptance-0.3.0.md`; no elevar resultados manuales ni prioridades sin evidencia identificable.
- [x] C.2 La matriz automatizada incluida en CI #39782991976 pasó: 417/417 tests relacionados con Clean2D; 19/19 amida/urea, 14/14 MOL P0, 12/12 importación estereoquímica SMILES, 13/13 Live MM stream, 5/5 proyección estereoquímica, 58/59 CompChem 3D, 35/36 geometry service y 71/76 casos ChemName acceptance. Los skips y límites quedan desglosados en `validation.md`; no equivale a aceptación manual universal.
- [x] C.3 Separar smoke automatizado, inspección estática y aceptación manual. La prueba Qt/OpenGL real, Windows interactivo, instalación y validación gráfica de los paquetes exactos se registran NOT TESTED.
- [ ] C.4 Ejecutar aceptación manual P0/P1 en los cuatro paquetes exactos de `a38fa804…` disponibles en Actions #39789295451, cuando el owner tenga acceso a Windows y CachyOS; registrar los resultados observados en la matriz existente sin inferir PASS de smoke automático.

## D. Preparación de publicación
- [x] D.1 Actualizar notas, OpenSpec y guía preview con SHAs, runs, artifact IDs, hashes y pasos de aceptación. La matriz `manual-acceptance-0.3.0.md` queda deliberadamente intacta; los resultados manuales pendientes no se modifican.
- [x] D.2 Dictamen: **candidata lista para aceptación manual; NO-GO para publicación**. Build Preview exacto `a38fa804…` está completo; aceptación manual P0/P1 y validación instalada/gráfica siguen pendientes. CI del código tiene riesgo de margen estrecho; detalles y evidencia en `validation.md`.
- [x] D.3 En el cierre original no se hizo merge, worktree adicional, tag, Release, PR ni cambio a canales públicos/`main`/`gh-pages`; ese estado se conserva como registro histórico. La autorización posterior del propietario se registra en el addendum de esta campaña.

## E. Integración final y preparación pública autorizadas (2026-10-23)
- [x] E.1 Confirmar fetch, árbol limpio, ancestro `61f7647b…`, cero divergencia y 29 commits de estabilización por delante; fast-forward-only desde `61f7647b…` hasta `d808070f…` y push normal de `release/v0.3.0-beta.1-prep`. `main` intacta.
- [x] E.2 Baseline proporcional, arquitectura, compileall, colecta, Ruff focal y 220 pruebas de seguridad/archivo/worker/release/version metadata; no ejecutar pytest masivo conforme a la política posterior.
- [x] E.3 Revisar flujos de archivo/proceso/AI/update, dependencias, permisos GitHub y matriz P0/P1; preparar cuerpo público curado y vincularlo en `release.yml` sin cambiar permisos, gates, canales ni lógica de publicación.
- [x] E.4 Build Preview #40144590463 de los cuatro paquetes desde `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`; verificar jobs, SHA común, provenance/checksums, exclusión Pillow, ChemName/RDKit y smoke real AppImage 2D→3D (detalles en validation.md).
- [ ] E.5 Completar los casos P0/P1 restantes en artifacts identificables; la aceptación limitada reportada el 2026-10-24 se documenta sin atribuirla a un SHA concreto ni promover casos no observados.
- [x] E.6 Configurar/verificar protección externa `v*` y entornos `beta`/`flatpak-beta` (completado 2026-10-24; evidencia y límites en `validation.md`). El tag se creó tras el Build Preview; el workflow Release quedó bloqueado por el skip marker documentado.

## F. Cierre público autorizado (2026-10-24)
- [x] F.1 Integrar por `git merge --ff-only` los cuatro commits de `release/v0.3.0-beta.1-security-prep`; documentar el alcance manual comunicado por el propietario sin adjudicar PASS a funciones no mencionadas.
- [x] F.2 Ejecutar auditoría OSV sobre la resolución actual Python 3.11 de los paquetes de build; cero advisories en 12 distribuciones. Dependabot/Code Scanning deshabilitados/no disponibles se informan como limitación, no como PASS.
- [x] F.3 Configurar ruleset activo para creación/actualización/borrado de `v*` y reviewer/ref policy `v*-beta.*` en `beta` y `flatpak-beta`; mantener permiso Actions predeterminado read.
- [x] F.4 Publicar la integración/documentación y verificar Build Preview #40144590463 para SHA `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`: cuatro paquetes, provenance/checksums, exclusión Pillow y gates funcionales automatizados.
- [x] F.5 Revisar tests focales de versión, workflow de release, sidecars AppImage/zsync, manifest/provenance y OpenSpec strict (**103 passed**, 2026-10-24); no ejecutar pytest completo.
- [ ] F.6 **BLOQUEADO:** el tag `v0.3.0-beta.1` apunta al SHA validado `85fa4ed2…`, pero el commit contiene `[skip ci]`; no apareció run de `release.yml` y el workflow no admite `workflow_dispatch` (API 422). No mover/borrar/recrear el tag; no hay Release/paquetes oficiales/Flatpak remoto. Requiere decisión del propietario sobre recuperación sin debilitar los gates.
- [ ] F.7 **BLOQUEADO:** resultados/hashes del preview y workflow ausente están registrados; no existen enlaces de Release ni assets oficiales para verificar. El commit documental de bloqueo deja la rama limpia sin alterar el tag SHA.
