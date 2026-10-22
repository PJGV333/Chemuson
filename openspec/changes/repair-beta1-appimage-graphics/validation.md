# Validación — reparar arranque gráfico AppImage beta.1

## Estado al inicio

- OpenSpec nuevo creado antes de tocar código: `repair-beta1-appimage-graphics`.
- Rama `fix/beta1-appimage-graphics` creada desde `stabilize/v0.3.0-beta.1` en `78b0139d95c5d5710e49c258f2355a5e95b7e478`. Un worktree principal; status limpio al inicio.
- No se movieron `main`, `release/v0.3.0-beta.1-prep` ni `release/v0.3.0-beta.1-rc-prep`; no se publicó Release.

## Baseline

- `compileall`: exit 0.
- Colecta errónea con `/usr/bin/python`: 2 errores de import RDKit (interpreter no apto); colecta correcta con `chem/bin/python`: 2598 tests, exit 0.
- Ruff global: un F401 basal (`math` sin usar en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`).
- Suite monolítica no repetida por el timeout previo de 8 minutos y la política de shards acotados del proyecto.

## Reproducción original — KDE/Wayland

- AppImage conservado en su ubicación original, SHA-256 `5559e28117d576b1b52249dc395aa42abd48b9cbc4bb2f84c483efdf8a643ea1`, tamaño 176645312 bytes.
- En la sesión real `XDG_SESSION_TYPE=wayland`, el AppImage se lanzó con `QT_DEBUG_PLUGINS=1`; a 8 s el proceso seguía vivo, pero no se observó ventana ChemUSON en screenshot del escritorio. El proceso se finalizó deliberadamente con SIGTERM, no ocurrió SIGABRT en esta repetición.
- Log: plugin `libqwayland`/integración Wayland EGL cargados desde `_MEIPASS`, después EGL no disponible, QOpenGLWidget no soportado y QRhi GLES2 sin contexto. Errores Fontconfig simultáneos contra `/etc/fonts/conf.d/48-guessfamily.conf`; causa aún sin atribuir.
- El log entregado por propietario acredita por separado el intento XWayland: la integración GLX del bundle no encuentra FBConfig y Qt no inicializa GLX. Se hará reproducción controlada XWayland.
- Evidencia exacta y SHA de log/captura en `baseline.md`; archivos guardados fuera del repo en `/tmp/chemuson-appimage-gl-repair-78b0139/`.

## Diagnóstico de causa raíz confirmado (antes de editar producción)

La causa de la pérdida de interfaz/OpenGL es el `libstdc++.so.6` incluido por PyInstaller, que es más antiguo que el runtime C++ requerido por Mesa del host. La ruta exacta es:

| Componente en la sesión CachyOS | AppImage original | `.venv` / host | Evidencia |
|---|---|---|---|
| Qt/plugins | PyQt6 metadata Qt 6.11.0 (`version=396032`); `libqwayland.so`, `libqt-plugin-wayland-egl.so`, `libqxcb.so` y `libqxcb-glx-integration.so` cargan desde `_MEIPASS` | PyQt6/Qt 6.11.2; plugins Qt de `chem/lib/python3.14/site-packages/PyQt6` | `QT_DEBUG_PLUGINS`, `LD_DEBUG`; los plugins existen y se cargan |
| GL/EGL/GLX/GLVND | `libGL.so.1.7.0`, `libEGL.so.1.1.0`, `libGLX.so.0.0.0`, `libGLdispatch.so.0.0.0`, `libGLX_mesa.so.0`, `libEGL_mesa.so.0`, `libgallium-26.2.4-arch3.1.so`, `libSPIRV-Tools.so`: `/usr/lib` del host | Mismos GL/EGL/GLX/Mesa del host | `/proc/<pid>/maps`, `LD_DEBUG=libs,versions`, `ldd`; drivers no están dentro del AppImage |
| Runtime C++ | `/tmp/_MEI…/libstdc++.so.6`, exporta máximo `GLIBCXX_3.4.30` | `/usr/lib/libstdc++.so.6.0.36`, exporta `GLIBCXX_3.4.36` | Mesa `/usr/lib/libSPIRV-Tools.so` requiere `GLIBCXX_3.4.32` |
| Fonts | `/tmp/_MEI…/libfontconfig.so.1` 2.13.1 y `libfreetype.so.6` 2.11.1 | `/usr/lib/libfontconfig.so.1.17.0` (Fontconfig 2.18.3), FreeType host 2.20.7 | `FcGetVersion()` sobre la copia exacta de `_MEIPASS`, `FT_Library_Version`, `fc-match --version`, `pkg-config`; mensajes de `/etc/fonts/conf.d/48-guessfamily.conf` |

Pruebas A/B de una variable:

1. Original, Wayland real, entorno heredado normal: plugins Wayland/EGL cargados, `EGL not available`, `QOpenGLWidget is not supported`, QRhi sin contexto y sin ventana visible. Original XCB/XWayland termina espontáneamente con código `-6` (SIGABRT); `LD_DEBUG` identifica en Mesa `GLIBCXX_3.4.32 not found` contra la copia PyInstaller con máximo `3.4.30`.
2. Quitar `LD_LIBRARY_PATH` mantiene el mismo error Wayland; la ruta heredada a Orca_QC no es la causa.
3. Pre-cargar sólo el `libstdc++.so.6` del host elimina la falla GL y la ventana ChemUSON aparece en screenshot XCB, pero los errores Fontconfig permanecen.
4. Pre-cargar sólo Fontconfig/FreeType del host suprime esos errores, pero mantiene el error `GLIBCXX_3.4.32` y el SIGABRT XCB. Por tanto los warnings Fontconfig son una segunda incompatibilidad, no la causa raíz GL.
5. Pre-cargar los tres ABI host produce ventana XCB visible y limpia el log de Fontconfig. No se usa como solución: valida qué debe ocurrir al dejar esos SONAMEs fuera del bundle.
6. El `.venv` muestra la ventana completa en KDE/Wayland y su sonda XCB obtiene contexto OpenGL 4.6 Core con renderer AMD/radeonsi; `LD_DEBUG` muestra `/usr/lib` para libstdc++, libfontconfig, freetype, GL/EGL/GLX. Esto confirma que no es un defecto del renderer ni de los shaders 3D.

Evidencia conservada bajo `/tmp/chemuson-appimage-gl-repair-78b0139/`: `original-xcb-ld/` (exit `-6` y trace fatal), `ab-clean-env/` (Wayland sin `LD_LIBRARY_PATH`), `ab-host-libstdc-xcb-window/` y `ab-host-desktop-libs/` (screenshots de ventana XCB), `ab-host-font-libs-only/` (Fontconfig arreglado pero GL sigue abortando), `bundled-font-libs/` (copias exactas del bundle/maps/versiones) y `venv-wayland/`, `venv-gl-reference/`. El AppImage original mantiene SHA-256 `5559e28117d576b1b52249dc395aa42abd48b9cbc4bb2f84c483efdf8a643ea1`.

## Corrección y verificación local (2026-10-22)

### Solución implementada

- `packaging/pyinstaller/runtime_filter.py` elimina de la tabla binaria de PyInstaller las SONAME `libstdc++.so.*`, `libfontconfig.so.*` y `libfreetype.so.*`. `chemuson.spec` ejecuta el filtro después de `Analysis` (incluidas dependencias transitivas) y antes de construir el EXE. No elimina Qt, RDKit, libgcc, Mesa, GL/EGL/GLX ni el renderer.
- `packaging/release/validate_appimage.py` rechaza esos SONAME en el CArchive y conserva los gates de iconos, RDKit, ChemName y startup offscreen; además exige una prueba gráfica XCB real con timeout de 180 s, contexto OpenGL Core ≥4.1, renderer llvmpipe, píxeles 2D y 3D, captura de ventana y cierre limpio.
- `src/chemuson/gui/packaged_graphics_smoke.py` usa la ventana/canvas reales, dibuja C-C-O con la herramienta de enlaces, abre el dock 3D, genera fixtures RDKit simple/complejo, comprueba framebuffer y cierra dock/ventana después del apagado de workers. La captura no es una aceptación manual.
- `build-preview.yml` instala `xvfb`, `xauth`, Mesa DRI, Fontconfig y FreeType para que el gate ejecute `xcb` sobre Xvfb/llvmpipe en el runner. Se registró el módulo GUI existente en `architecture/modules.yml`. Se añadieron regresiones del filtro/validador y una corrección al inventario ChemName: el CArchive contiene `templates/__init__.py` como dato adicional, por lo que el inventario contractual cuenta sólo recursos `.mol` y no falsamente ese archivo Python.

### Pruebas realizadas

- Smoke fuente con `chem/bin/python`, `DISPLAY=:0`, `QT_QPA_PLATFORM=xcb`, `LIBGL_ALWAYS_SOFTWARE=1`, `GALLIUM_DRIVER=llvmpipe`: PASS. Ventana visible, canvas 2D (3 átomos/2 enlaces, 684 muestras modificadas), contexto OpenGL 4.6 Core/llvmpipe, escena simple (9 átomos, 11.376 vértices, 4.703 muestras) y compleja (24 átomos, 31.464 vértices, 4.862 muestras), dock/ventana cerrados sin SIGABRT. Captura: `/tmp/chemuson-appimage-gl-repair-78b0139/source-graphics-smoke-v9/window.png` (95 KB), verificada visualmente.
- PyInstaller local 6.22.3 produjo un EXE de prueba de 321 MB con `chem/bin/python` 3.14. `_validate_pyinstaller_resources` inspeccionó el CArchive real y pasó módulos Qt, 69 SVG, nueve MOL y ausencia de los tres ABI host. El reporte `graphics-smoke.json` del binario congelado pasó y la captura muestra 2D y el modelo 3D complejo: `/tmp/chemuson-appimage-gl-repair-78b0139/frozen-build-py314/validation/graphics-smoke.png` (97 KB). Reporte/logs en el mismo directorio.
- La invocación local de `_validate_frozen_graphics` usó un wrapper `xvfb-run` temporal que mantuvo el `DISPLAY=:0` existente. Esto valida CLI congelado, CArchive, XCB, software Mesa y lógica del validador, pero **no** equivale a Xvfb real. En este host no hay `xvfb-run`, `Xvfb` ni `appimagetool`; por ello no se construyó un AppImage nuevo ni se probó el comando oficial end-to-end. Build Preview instalará Xvfb en CI.
- Verificación final agrupada: `chem/bin/python -m pytest -q tests/test_appimage_validation.py tests/test_linux_distribution_scripts.py tests/test_viewport3d_widget.py tests/architecture/test_import_boundaries.py tests/architecture/test_module_catalog.py tests/test_ci_pytest_shards.py tests/test_test_workflow_contract.py`: **159 passed**.
- Se actualizó `tests/ci/expected_pytest_nodeids.txt` para los seis tests nuevos. El plan real de colección (`tools/ci_pytest_shards.py plan --shards 8 ...`) coincide exactamente con **2.604 IDs**, ocho shards, estimación máxima 208,144 s; el hash del plan es `f71bbf819fcbfb90066f3bf0a9f61a6006e5eff605f2df9bf10e4693f94a6030`.
- Actions [#39809759478](https://github.com/PJGV333/Chemuson/actions/runs/39809759478), exact SHA `db48d71574850dc7c7c72961e34f5fc782ff3e4a`: plan 2.604 IDs y los ocho shards terminaron con exit 0, sin timeout y con JUnit/timings para todos: **2.582 passed, 22 skipped, 0 fallos/errores**. El shard más lento fue 5 (295,104 s de 300). Windows smoke pasó. El run global quedó **cancelled**, no success: `flatpak-smoke` agotó sus 10 min dentro de `apt-get install flatpak-builder`; el log muestra descargas Ubuntu extremadamente lentas y el job se canceló antes de la validación Flatpak. No es un fallo de los tests ni del código Flatpak modificado (no se cambió); no se repitió el job exacto.
- `chem/bin/python -m compileall -q src tests tools packaging`, Ruff focal (`F401,F811,F821,E722,E741` en archivos Python afectados), `openspec validate repair-beta1-appimage-graphics --strict` y `git diff --check`: PASS. Ruff global devolvió únicamente la excepción basal `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se tocó.
- No se repitió la suite completa local ni ningún shard/job que hubiese agotado su límite; el run remoto ya cubrió 2.604 IDs en shards limitados a 300 s. Se usa además el baseline #39782991976 previamente documentado.

### Límites de entrega y aceptación

- El EXE PyInstaller es una verificación de paquete local, construido con Python 3.14 en CachyOS; no representa el build oficial Ubuntu 22.04/Python 3.11 ni un artefacto AppImage con provenance. No se ha creado Build Preview para el SHA corregido.
- No se repitió la prueba del candidato exacto en KDE/Wayland, no se generó el AppImage corregido y no hay Xvfb real local. La prueba software XCB/llvmpipe es automatizada y no sustituye la instalación ni aceptación del propietario.
- Se preservaron el AppImage original (SHA-256 `5559e28117d576b1b52249dc395aa42abd48b9cbc4bb2f84c483efdf8a643ea1`) y las referencias `main`, `release/v0.3.0-beta.1-prep` y `release/v0.3.0-beta.1-rc-prep`; no se hizo merge, PR, tag, Release ni publicación. `docs/release/manual-acceptance-0.3.0.md` permanece intacto; ítems no probados siguen `NOT TESTED`.

## Dictamen

**NO-GO para publicación.** La corrección de ABI y el smoke gráfico source/frozen en XCB+llvmpipe pasan localmente; quedan pendientes el AppImage exacto con provenance, el gate ejecutado realmente bajo Xvfb/runner, la repetición Wayland del candidato y la aceptación manual P0/P1.

## Build Preview exacto — validación final 2026-10-22

- Se confirmó `origin/fix/beta1-appimage-graphics` en `fa51c896992c2d72d01ebf03deca316d2cce175d` y worktree limpio antes de iniciar. Se creó y publicó `release/v0.3.0-beta.1-graphics-rc-prep` apuntando al mismo SHA, sin commit nuevo. Las referencias de release preexistentes quedaron intactas: `release/v0.3.0-beta.1-prep`=`61f7647b683db3eca19ea06e5ec82ab2897d640a`; `release/v0.3.0-beta.1-rc-prep`=`a38fa8045c858f8ff4b6002950f57e50292aea81`.
- Build Preview oficial [#39848401263](https://github.com/PJGV333/Chemuson/actions/runs/39848401263), evento `workflow_dispatch`, validó el branch y SHA exactos. Identidad: versión `0.3.0-beta.1`, rama `release/v0.3.0-beta.1-graphics-rc-prep`, SHA fuente `fa51c896992c2d72d01ebf03deca316d2cce175d`, `publication=false`. El reporte descargado declara `complete=false`, `build_status=failure`: source validation PASS; Windows portable/installer PASS; Flatpak PASS; AppImage FAIL.
- El runner sí generó el Type 2 llamado `Chemuson-v0.3.0-beta.1-preview-fa51c896-linux-x86_64.AppImage`. El digest SHA-256 que imprimió su verificador fue `47042dbc20e59fe6f9b4bb5753796b4552637b5fc7178d0f2e45d1b79053c28a`; es evidencia de un archivo efímero del runner, no un artefacto descargable/aceptable. El job falló antes del manifiesto `preview-provenance.json`, checksums y upload; **no existe artifact AppImage ni enlace de descarga** para este SHA.
- Pasó la inspección CArchive, incluido el rechazo de `libstdc++.so.*`, `libfontconfig.so.*`, `libfreetype.so.*`; el ejecutable congelado también pasó iconos (69 SVG, tema claro/oscuro/HiDPI), RDKit/worker 3D, ChemName (9 templates y 6 moléculas), y el plugin `xcb` fue encontrado. El gate de GUI real en `xvfb-run` falló al inicializar el plugin Qt `xcb`, con mensaje “xcb-cursor0 or libxcb-cursor0 is needed”; el mismo job muestra `libxcb-cursor0` instalado correctamente por apt. La causa exacta de resolución/carga dinámica queda **sin determinar**. Por eso el test gráfico no pasa y no se afirma que se hayan validado todas las dependencias requeridas ni la pantalla/cierre.
- Artifacts existentes del run (ninguno es AppImage): Windows portable [12455514065](https://github.com/PJGV333/Chemuson/actions/runs/39848401263/artifacts/12455514065), Windows setup [12457180560](https://github.com/PJGV333/Chemuson/actions/runs/39848401263/artifacts/12457180560), Flatpak [12459202406](https://github.com/PJGV333/Chemuson/actions/runs/39848401263/artifacts/12459202406), reporte [12459311217](https://github.com/PJGV333/Chemuson/actions/runs/39848401263/artifacts/12459311217). Sólo el reporte acredita la identidad/estado del workflow; no sustituye provenance del AppImage inexistente.

### Investigación separada — cancelación de Flatpak

- El job `flatpak-smoke` del workflow `test` [#39809759478](https://github.com/PJGV333/Chemuson/actions/runs/39809759478), SHA `db48d71574850dc7c7c72961e34f5fc782ff3e4a` (no el SHA fuente del preview actual), terminó `cancelled` a los 10 minutos mientras seguía `sudo apt-get install -y flatpak-builder`; no alcanzó la validación del manifest. Apt enumeró 123 paquetes (99,0 MB), con descargas largas/intermitentes desde `azure.archive.ubuntu.com`, entradas `Ign:` y cancelación cuando iba por los paquetes 46/123; el `apt-get update` había completado. Evidencia/log completo: `/tmp/chemuson-beta1-graphics-final-baseline/flatpak-cancelled-test-job.log` (SHA-256 `cdffe2928e18284ba4e67ae7642da49fe056fe15d9fd54b0cf4ad258b09704ea0`).
- En el propio Build Preview #39848401263, la instalación de dependencias Flatpak pidió de nuevo 123 paquetes/99 MB y tardó aproximadamente 38 min 23 s; después el bundle se compiló, validó y subió correctamente. El job completo duró 43m06s (límite existente 45m). Esto respalda una limitación/intermitencia de infraestructura de paquetes, no un fallo funcional Flatpak; **no se aumentó ningún timeout**. La cancelación del test workflow y el PASS del job Flatpak Build Preview son resultados distintos.

### Estado de aceptación y evidencia preservada

- No se abrió el AppImage anterior ni se repitieron diagnósticos manuales del paquete viejo. No se probaron manualmente este candidate ni UI en CachyOS/KDE/Wayland. Inicio, editor 2D, Clean2D, visor 3D, cierre limpio, XWayland interactivo y aceptación manual permanecen `NOT TESTED`.
- Logs completos del AppImage, Flatpak preview, Flatpak cancelado, reporte Build Preview y baseline local se guardaron bajo `/tmp/chemuson-beta1-graphics-final-baseline/`; los hashes de job logs están en esa carpeta y en el informe de agente.
- No hubo merge, tag, GitHub Release ni publicación. La matriz manual sigue sin modificaciones de estado; los casos sin evidencia de aceptación del propietario permanecen `NOT TESTED`.

## Seguimiento: carrera del smoke y reparación del harness — 2026-10-22

- Build Preview [#39867831993](https://github.com/PJGV333/Chemuson/actions/runs/39867831993), source SHA `7b1c450321567fbe4903393ee530be87caee102c`, validó `publication=false`. Source, Windows portable/setup y Flatpak pasaron. AppImage Type 2 se generó, pero el gate gráfico falló en la escena compleja: `scene_atoms=24`, `scene_bonds=25`, contexto GL válido (4.5 Core, llvmpipe, Mesa 23.2.1), pero `scene_empty`, renderer CPU con 0 vértices y `changed_samples=410`; falló antes de provenance/checksum/upload. No existe AppImage artifact ni enlace de descarga.
- Causa demostrada: al dibujar la fixture 2D, `MainWindow` programa la sincronización 3D con debounce de 180 ms. El trabajo automático arrancaba mientras el smoke inyectaba/renderizaba escenas 3D propias y vaciaba la escena de viewport. El fallo cambia de escena simple a compleja según el timing; no es fallo de inicialización XCB, GL ni creación del mesh. Se corrigió únicamente el harness: detiene el timer y suprime el scheduling de geometría automática en la instancia efímera del smoke. No se cambia el comportamiento de la aplicación ni de la geometría química. Se añadió regresión focal y un node ID esperado.
- Reproducción tras la corrección: EXE congelado local con Python 3.11/CI dependencies, Xvfb Ubuntu 21.1.4 y Mesa 23.2.1/llvmpipe alcanzó OpenGL 4.5 Core; escena simple 11.376 vértices/2.849 muestras modificadas y compleja 31.464/2.794; canvas 2D visible y dock/ventana cerrados limpiamente. `_validate_frozen_graphics` también pasó con Xvfb Ubuntu y Mesa 26 del host, sin heredar `LD_LIBRARY_PATH`. Artefactos de prueba local: `/tmp/chemuson-ci-compatible-py311/fixed-ubuntu-mesa-report.json`, `fixed-ubuntu-mesa.png`, `fixed-ubuntu-mesa-run.log`; no son AppImage de preview.
- Test workflow [#39867829904](https://github.com/PJGV333/Chemuson/actions/runs/39867829904), mismo SHA diagnóstico `7b1c450…`: **2.606/2.606 node IDs**, 2.584 passed, 22 skipped, 0 fallos/errores; ocho shards completos y máximo 262,03 s, Windows smoke PASS. Ese run precede el nuevo test del harness. Con el test añadido, colección local = 2.607 IDs, plan 8 shards válido, hash `02cb5ad4442c263c6218b0517cde9c82db4c9a46f49149bfa52f19d992d4d305`, máximo estimado 208,219 s.
- Verificación local de cambios: packaging/Linux focus **25 passed**, arquitectura **283 passed**, workflow/sharding contracts **10 passed**; `compileall`, Ruff focal, OpenSpec strict y `git diff --check` PASS. Ruff global conserva sólo el F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. El test nuevo y la colección 2.607 deben confirmarse en el CI del SHA corregido.
- Estado: se necesita un Build Preview nuevo sobre el harness corregido. Hasta que valide y suba provenance/checksums, no hay artifact AppImage para CachyOS. No se ha realizado aceptación manual; inicio, editor 2D, Clean2D, visor 3D y cierre siguen `NOT TESTED`. No se alteraron timeouts ni se publicó nada.
