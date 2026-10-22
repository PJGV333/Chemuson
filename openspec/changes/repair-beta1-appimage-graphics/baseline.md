# Baseline — reparación gráfica AppImage beta.1

Fecha: 2026-10-22. Repositorio `/home/ccachyavgp/Documentos/ChemUSON`.

## Git y preservación

- Rama inicial `stabilize/v0.3.0-beta.1`, HEAD `78b0139d95c5d5710e49c258f2355a5e95b7e478` (`[skip ci] docs(beta): record RC preview audit`); árbol limpio antes de abrir este OpenSpec. La reparación parte de ese HEAD y preserva el único worktree principal.
- Referencia de build defectuosa no se mueve: `origin/release/v0.3.0-beta.1-rc-prep` permanece `a38fa8045c858f8ff4b6002950f57e50292aea81`; `main` `f0fde603371255902bf0e630ca1ce032b8f16ad3` y `release/v0.3.0-beta.1-prep` `61f7647b683db3eca19ea06e5ec82ab2897d640a` no se cambian.
- Se creó la rama local `fix/beta1-appimage-graphics` desde `78b0139d…`; ningún otro worktree. El AppImage original se deja en la ruta de Descargas sin escribirlo ni sustituirlo.
- AppImage original `/mnt/HDD_4TB/Descargas_provisional/chemuson-preview-linux-appimage/Chemuson-v0.3.0-beta.1-preview-a38fa804-linux-x86_64.AppImage`: 176,645,312 bytes, SHA-256 `5559e28117d576b1b52249dc395aa42abd48b9cbc4bb2f84c483efdf8a643ea1`; coincide con el artifact AppImage de Actions #39789295451.
- Log original aportado `/tmp/chemuson-gl-debug.log`: 24,400 bytes. Se conserva intacto. Evidencias nuevas en `/tmp/chemuson-appimage-gl-repair-78b0139/`.

## Estado inicial de tests/linters

- `git status --short`: vacío; estado exacto `/tmp/chemuson-appimage-gl-repair-78b0139-baseline/git-status-short.txt`.
- `python -m compileall src tests tools packaging`: exit 0; salida en `compileall.txt`.
- `pytest --collect-only -q` con `/usr/bin/python` (Python 3.14.8 de sistema): exit 2, `2572 tests collected, 2 errors`; las dos colecciones fallidas importan `rdkit` (`tests/test_mol_interoperability_p0.py`, `tests/test_smiles_stereo_import.py`). Es un intérprete incorrecto para el proyecto, no se atribuye al código.
- `chem/bin/python -m pytest --collect-only -q`: exit 0, `2598 tests collected in 0.68s`; Python 3.14.8, RDKit 2026.03.6 y PyQt6 disponibles desde el venv del repo. Salidas de ambos intérpretes conservadas en el directorio baseline.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 sólo por F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no relacionado.
- No se repitió `pytest -q` monolítico: el OpenSpec de estabilización conserva su intento de 8 min sin resumen, y la skill/protocolo bounded ordena no repetir el agregado tras timeout. Los focos gráficos y de empaquetado se ejecutarán con timeout externo <10 min.
- Baseline raw completo: `/tmp/chemuson-appimage-gl-repair-78b0139-baseline/` (incluye salida del intérprete sistema y la colecta válida del repo). SHA de los principales logs se registrará al cierre.

## Entorno real y reproducción del paquete defectuoso

- Sesión disponible: `DISPLAY=:0`, `WAYLAND_DISPLAY=wayland-1`, `XDG_SESSION_TYPE=wayland`, `XDG_CURRENT_DESKTOP=KDE`; KWin Wayland, Plasma y XWayland activos. `/dev/dri` y `glxinfo`/`eglinfo` están disponibles. `glxinfo -B` del anfitrión devuelve AMD Ryzen 9 7950X/radeonsi, Mesa 26.2.4, OpenGL 4.6 y `direct rendering: Yes`.
- Se lanzó el AppImage original en la sesión Wayland con el entorno gráfico normal y `QT_DEBUG_PLUGINS=1`. Tras 8 s el proceso seguía vivo pero no había ventana ChemUSON visible en la captura `reproduction/original-wayland.png`; la prueba lo terminó deliberadamente con SIGTERM (exit 143), no se contabiliza como SIGABRT espontáneo.
- El log de esa reproducción muestra metadatos Qt del plugin `libqwayland.so` y que carga desde el bundle `libqt-plugin-wayland-egl.so`, pero después: `qt.qpa.wayland: EGL not available`, `QOpenGLWidget is not supported on this platform`, `QRhiGles2: Failed to create temporary context`, `QRhiGles2: Failed to create context`, `Failed to create QRhi for QBackingStoreRhiSupport`. Evidencia: `reproduction/original-wayland.log` y `original-wayland.status`.
- El mismo log registra errores Fontconfig en `/etc/fonts/conf.d/48-guessfamily.conf` (`invalid attribute 'xsi:nil'`) y warnings `invalid constant used`; también warnings en `48-spacing.conf` y `49-sansserif.conf`. Esto sucede durante el mismo startup pero aún no está demostrado que cause la ausencia de ventana.
- El log histórico del propietario termina con `QXcbGlxIntegrationPlugin` del bundle cargado y tres fallos `qglx_findConfig: Failed to finding matching FBConfig`, seguidos por `Could not initialize GLX`; es evidencia XWayland separada, aún debe reproducirse de forma controlada.
- La captura inicial muestra el navegador de fondo, no ChemUSON. Por ello el resultado no se etiqueta como PASS visual ni como fallo de cierre de producto; es una reproducción de startup Wayland sin ventana visible.

## Actualización diagnóstica posterior al baseline

Sin cambiar el AppImage original ni el host, se completó la inspección dinámica antes de editar producción. El detalle de rutas, versiones, A/B y sus límites está en `validation.md` y `design.md`. Conclusión confirmada: conflicto entre Mesa del host (requiere `GLIBCXX_3.4.32`) y `libstdc++` congelado (máximo `3.4.30`); Fontconfig/FreeType congelados (2.13.1/2.11.1) generan una segunda incompatibilidad de sintaxis con configuración CachyOS moderna. La causa principal fue reproducida en Wayland/XWayland y atribuida con preload A/B; se actualiza el plan de reparación antes de cambiar código.

El entorno contiene `LD_LIBRARY_PATH=/home/ccachyavgp/Aplicaciones/Orca_QC` duplicado, pero retirarlo no cambia el fallo. La reproducción original XCB exit `-6` es espontánea y separada de las pruebas donde la terminación SIGTERM se hizo deliberadamente. Los logs/screenshot/maps/copias exactas del runtime están en `/tmp/chemuson-appimage-gl-repair-78b0139/`; su SHA-256 se recopilará al cerrar el cambio.

## Baseline de validación final del paquete (2026-10-22)

- HEAD local y `origin/fix/beta1-appimage-graphics`: `fa51c896992c2d72d01ebf03deca316d2cce175d`; worktree limpio antes de iniciar. Se confirmó además que el remoto `fix/beta1-appimage-graphics` apunta exactamente a ese SHA.
- `python -m compileall src tests tools packaging`: exit 0.
- `pytest --collect-only -q`: exit 0, 2.604 tests collected in 0.79s.
- `pytest -q`: el comando excedió el límite externo de 720s; la última salida mostró progreso de 41%, sin resumen final. No se atribuye PASS ni fallo de tests; no se repitió monolíticamente.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por el F401 basal `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; sin otros hallazgos.
- Salida exacta capturada en `/tmp/chemuson-beta1-graphics-final-baseline/exact-output.txt` (SHA-256 `8184301ffb1a33f6c6eb378f4c152afe3d7f869b548a94fd1e3cbb77913909ad`). La captura incluye status, compileall, colección, progreso parcial del agregado y salida Ruff; el timeout lo informó la herramienta a los 720s.
