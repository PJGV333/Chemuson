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

## Corrección / validación

Pendiente. Actualizar con hipótesis probadas, cambios, tests, re-build exacto, capturas GUI, resultado de visor 3D y límites manuales antes de cerrar.

## Dictamen provisional

**NO-GO para publicación.** La reproducción Wayland no muestra ventana. No se ha modificado producción ni se ha generado candidato corregido. La rama RC defectuosa permanece intacta.
