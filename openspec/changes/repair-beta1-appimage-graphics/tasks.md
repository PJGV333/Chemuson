# Tareas — reparar el AppImage beta.1

## 1. Baseline y OpenSpec
- [x] 1.1 Capturar estado/HEAD/refs, compileall, colección del intérprete `chem/bin`, Ruff focal y límites de la suite; registrar la incompatibilidad del `pytest` de sistema que carece de RDKit.
- [x] 1.2 Preservar SHA/tamaño del AppImage original y copia del log entregado; reproducirlo sin cambiar sus bytes.
- [x] 1.3 Completar causa/procedencia/historial en `validation.md`; actualizar `design.md` con evidencia antes de tocar producción.

## 2. Reproducción y diagnóstico dinámico
- [x] 2.1 Reproducir startup original en KDE/Wayland real y XWayland con logs, estado de PID y captura; distinguir salida involuntaria de terminación de prueba.
- [x] 2.2 Extraer AppImage a `/tmp`, revisar `AppRun`, plugins QPA, dependencias ELF/SONAME/versiones, PyInstaller archive y rutas temporales de bibliotecas.
- [x] 2.3 Capturar `LD_DEBUG=libs,versions`, `QT_DEBUG_PLUGINS`, `/proc/<pid>/maps`, entorno, `readelf -d/-V` y `ldd` sobre binario/plugins/bibliotecas sospechosas; comparar la app `.venv` en el mismo host.
- [x] 2.4 Demostrar o descartar separadamente conflicto libGL/libEGL/libGLX/GLVND/Mesa/DRM, dependencia ausente, procedencia/incompatibilidad de plugins, inicialización prematura del visor/QRhi y errores Fontconfig/FreeType.
- [x] 2.5 Documentar la causa raíz y una prueba A/B controlada que cambie una sola condición; si la causa no se demuestra, no modificar producción por conjetura.

## 3. Reparación mínima
- [x] 3.1 Implementar sólo la solución respaldada por las trazas. No cambiar el stack del host, los archivos de configuración del usuario, los drivers ni deshabilitar OpenGL/3D.
- [x] 3.2 Mantener el editor 2D independiente de inicialización 3D no solicitada, sólo si se demuestra que la creación del visor bloquea la ventana; la comparación A/B descartó esa hipótesis, por lo que no se cambia el ciclo de vida del visor.
- [x] 3.3 Investigar/mitigar Fontconfig independientemente; excluir el par Fontconfig/FreeType bundled incompatible y no silenciar advertencias ni reemplazar config del sistema.
- [x] 3.4 Añadir regresiones focalizadas de filtrado ABI del host y validación fail-closed del archive, sin alterar química ni serialización.

## 4. Validación de paquete y GUI
- [x] 4.1 Añadir al validador AppImage un startup de plataforma gráfica real con xcb+Xvfb+Mesa software/llvmpipe (no `offscreen`), captura de ventana visible, timeout y comprobación de cierre/no SIGABRT; mantener gates offscreen previos. Build Preview #398678 confirmó inicialización XCB y contexto OpenGL 4.5 Core/llvmpipe bajo Mesa 23.2.1; identificó una carrera del harness, separada del runtime del paquete, ahora suprimida sólo en la instancia de smoke.
- [x] 4.2 Añadir cobertura de renderer que abra/cierre el dock 3D y renderice moléculas simple y compleja con OpenGL/Mesa; demostrar iluminación y renderizado no degradados en el smoke software local.
- [ ] 4.3 Pendiente de rerun oficial sobre el harness corregido: Build Preview #39867831993, SHA `7b1c450321567fbe4903393ee530be87caee102c`, generó un AppImage Type 2 efímero y el contexto XCB/OpenGL fue válido, pero el render de la escena compleja quedó vacío antes de provenance/checksum/upload. Diagnóstico: el timer de sincronización 3D de la fixture 2D reemplazaba asíncronamente las escenas inyectadas. El harness ahora detiene el timer y bloquea el scheduling; la prueba congelada local con Xvfb/Mesa 23.2.1 pasa, pero aún no hay artifact AppImage descargable.
- [ ] 4.4 Reprobar original/candidato corregido en Wayland KDE real, XWayland y modo Mesa software, con evidencia visible antes/después y cierre limpio.
- [x] 4.5 Ejecutar pruebas focalizadas, arquitectura si aplica, OpenSpec estricto, compileall, Ruff y diff-check. No repetir suite/shard exacto que excedió límites; documentar excepciones. En el SHA diagnóstico `7b1c450`, Actions #39867829904 verificó 2.606/2.606 node IDs: 2.584 passed, 22 skipped, 0 fallos; el test nuevo del harness se añade a la colección 2.607 y requiere su propio run CI.

## 5. Entrega
- [x] 5.1 Actualizar OpenSpec/AGENT_REPORT con archivos, causa, bibliotecas/procedencia, pruebas y limitaciones.
- [x] 5.2 Ejecutar Build Preview oficial desde el ref `release/v0.3.0-beta.1-graphics-rc-prep`; #398484 (`fa51c896…`) falló al cargar XCB; #39867831993 (`7b1c450…`, `publication=false`) ya inicializó XCB pero falló por la carrera del harness. Ninguno produjo artifact AppImage.
- [ ] 5.3 Commit pequeño y push normal únicamente a `fix/beta1-appimage-graphics`, actualizar el ref RC autorizado al mismo SHA y verificar remoto/árbol limpio; después ejecutar Build Preview y workflow de tests sobre ese SHA.
- [x] 5.4 Sin merge, PR, tag, GitHub Release o publicación. NO-GO se mantiene hasta aceptación del propietario en el AppImage exacto.
