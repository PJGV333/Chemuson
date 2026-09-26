# Diagnóstico KDE real

Capturado antes de modificar código, desde el `.venv` de este checkout, sin forzar `QT_QPA_PLATFORM`.

- `XDG_SESSION_TYPE=wayland`
- `WAYLAND_DISPLAY=wayland-0`
- `DISPLAY=:0` (XWayland está disponible)
- `QT_QPA_PLATFORM=<unset>`
- `QT_SCALE_FACTOR=<unset>`
- `QT_AUTO_SCREEN_SCALE_FACTOR=<unset>`
- `QT_SCREEN_SCALE_FACTORS=<unset>`
- Plataforma Qt: `wayland`
- Pantalla primaria: `HDMI-A-4`
- DPR: `2.0`
- DPI lógico: `96.0`
- DPI físico: `40.80770812061806`
- Geometría: `(0, 0, 1422, 800)`
- Geometría disponible: `(0, 0, 1422, 800)`

Python activo: `/home/ccachyavgp/Documentos/ChemUSON-UI/.venv/bin/python`
Módulo importado: `/home/ccachyavgp/Documentos/ChemUSON-UI/src/chemuson/__init__.py`

## Reproducción visual requerida

La sesión usa Wayland, DPR 2.0. Las ejecuciones normales de la aplicación y del spike, y la captura directa con Spectacle, deben registrarse después de abrir las ventanas. No se acepta `offscreen` como verificación final.
