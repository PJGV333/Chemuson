# Observaciones de cierre — modernización UI

Este registro conserva únicamente asuntos vigentes o límites verificables. El contexto y la historia de los ya resueltos están en [`docs/history/CAMPAIGNS.md`](../history/CAMPAIGNS.md); los fallos no se ocultan ni se convierten en excepciones.

## Baseline de QA (no atribuible a la UI)

- Suite completa al cierre: 1760 passed, 55 skipped y un fallo preexistente en `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`.
- Ruff scoped mantiene el F401 de `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; Clean2D quedó fuera del alcance UI.
- El plugin Qt `offscreen` imprime `This plugin does not support propagateSizeHints()` durante algunos smokes. Es un warning no fatal; no se observó excepción funcional.

Estos puntos requieren alcance separado; no cambiar su baseline ni su implementación desde esta ficha.

## Límites de alcance

- La aceptación visual KDE/Wayland de la UI integrada se registró como aprobada.
- Clean2D y la química/geometría de plantillas (Fischer, Haworth, conformaciones de silla, tetrandrina y otros casos delicados) no se alteraron en las fases UI.
- No se produjo montaje before/after: no hay captura original homogénea de Fase 0 versionada.
- El flujo portable existente llamado `.AppImage` no fue validado como contenedor AppImage Type 2: faltaba `appimagetool`. Flatpak y el arranque PyInstaller sí pasaron los smokes descritos en el cierre.

No hay un issue visual UI abierto registrado al cierre. Cualquier nueva observación necesita reproducción, clasificación respecto a baseline y un OpenSpec propio.
