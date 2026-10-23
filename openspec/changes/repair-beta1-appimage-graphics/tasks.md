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
- [x] 4.1 Añadir al validador AppImage un startup de plataforma gráfica real con xcb+Xvfb+Mesa software/llvmpipe (no `offscreen`), captura de ventana visible, timeout y comprobación de cierre/no SIGABRT; mantener gates offscreen previos. Build Preview #39869114647 pasa el gate con XCB, OpenGL 4.5 Core/llvmpipe/Mesa 23.2.1, canvas 2D visible, dos escenas 3D y cierre limpio. El harness suprime su timer automático para no reemplazar las escenas de prueba.
- [x] 4.2 Añadir cobertura de renderer que abra/cierre el dock 3D y renderice moléculas simple y compleja con OpenGL/Mesa; demostrar iluminación y renderizado no degradados en el smoke software local.
- [x] 4.3 Validar Build Preview #39869114647, source SHA `7c7a02c51c2af5d32b320839b77d544a0271af1c`: workflow completo, `publication=false`, artifact AppImage Type 2 [12465108715](https://github.com/PJGV333/Chemuson/actions/runs/39869114647/artifacts/12465108715). Nombre `Chemuson-v0.3.0-beta.1-preview-7c7a02c5-linux-x86_64.AppImage`, 175.350.976 bytes, SHA-256 `6f9e8205820c8662fca316fb9312590d2236ab29d6ba93a5ad6689922c26619e`; checksums y provenance verificados. Gate real Xvfb/XCB/Mesa pasó con ambas escenas; artifact descargado y extraído sin FUSE, `AppRun --version` y smoke gráfico congelado local también pasan.
- [ ] 4.4 Reprobar original/candidato corregido en Wayland KDE real, XWayland y modo Mesa software, con evidencia visible antes/después y cierre limpio.
- [x] 4.5 Ejecutar pruebas focalizadas, arquitectura si aplica, OpenSpec estricto, compileall, Ruff y diff-check. Actions test #39869114633, mismo SHA candidato, verificó 2.607/2.607 node IDs: 2.585 passed, 22 skipped, 0 fallos; ocho shards completos (máximo 269,215 s), Windows y Flatpak smoke PASS. Un run duplicado concurrente #39869114496 perdió el shard 5 por el límite de 300 s (422/439 ejecutados) y fue cancelado en apt Flatpak; el mismo SHA tiene el run completo verde y no se repitió el shard fallido.

## 5. Entrega
- [x] 5.1 Actualizar OpenSpec/AGENT_REPORT con archivos, causa, bibliotecas/procedencia, pruebas y limitaciones.
- [x] 5.2 Build Preview oficial #39869114647 desde `release/v0.3.0-beta.1-graphics-rc-prep`, source SHA exacto `7c7a02c51c2af5d32b320839b77d544a0271af1c`; reporte completo/success, AppImage, Windows y Flatpak PASS, `publication=false`. Los intentos #398484 y #398678 permanecen como evidencia histórica de fallos ya corregidos.
- [x] 5.3 Commit `7c7a02c51c2af5d32b320839b77d544a0271af1c` y push normal sólo a `fix/beta1-appimage-graphics` y al ref RC creado/permitido `release/v0.3.0-beta.1-graphics-rc-prep`; ambos apuntaron al SHA exacto del paquete. Refs de release preexistentes intactas.
- [x] 5.4 Sin merge, PR, tag, GitHub Release o publicación. NO-GO se mantiene hasta aceptación del propietario en el AppImage exacto.

## 6. Causa raíz del visor 3D vacío en el AppImage (flujo normal 2D→3D)

- [x] 6.1 Reproducir el mismo escenario en fuente y binario congelado (dibujar 2D → solicitar 3D por el flujo normal del servicio) y clasificar la falla: R1 fuente ok=true; R2 worker de flujo continuo congelado exit 2 por argv inválidos; R3 modo interno no-stream ya funcionaba (RDKit empaquetado sano); R4 eventos NDJSON iban a stdout en vez del archivo de respuesta. Causa: transporte de flujo continuo congelado roto + enmascaramiento de estado en `paintGL`.
- [x] 6.2 Corregir solo el transporte/empaquetado: `_EVENT_STREAM` en `_rdkit_worker.py`, rama `_run_worker_stream_frozen` en `rdkit_safe.py`, guardado de `render_requested` a escenas no vacías en `viewport3d/widget.py`. Sin cambios químicos ni de serialización.
- [x] 6.3 Smoke de flujo normal congelado end-to-end (`geometry3d/flow_smoke.py` + flag oculto), catalogado en OpenSpec/arquitectura; validadores de paquete lo ejecutan y exigen `ok=true`.
- [x] 6.4 Variante de build de diagnóstico separada de la distribución pública: `packaging/release/diagnostic_build.py` + `BUILD_TYPE=diagnostic`; conserva integridad, provenance y SHA-256; nunca publica ni actualiza; nombre marcado; validador fail-closed lo verifica.
- [x] 6.5 Pruebas focalizadas de regresión: 3 pruebas nuevas de worker congelado en `tests/test_rdkit_live_stream.py` (relay ordenado de eventos, salida distinta de cero sin eventos, enrute al flujo de salida explícito); suites focales y arquitectura en verde.
- [ ] 6.6 Rebuild PyInstaller limpio sin artefactos ni cachés previos y validadores empaquetados (iconos, RDKit worker, ChemName, flujo geométrico) sobre el binario nuevo.
- [ ] 6.7 AppImage de diagnóstico real Type 2 en CachyOS vía `build_appimage.sh` `BUILD_TYPE=diagnostic`, validación completa con geometría 3D del flujo normal y manifiesto de provenance + SHA-256.
- [ ] 6.8 Aceptación funcional separada del build: flujo normal 2D→conformero→3D desde el AppImage extraído y desde Dolphin en KDE; aceptación visual final del propietario.
