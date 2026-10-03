# Modernización de UI de ChemUSON

**Estado:** Fases 0–8 cerradas. Fase 9 no iniciada.<br>
**Versión:** se conserva `0.3.0-dev`; no hubo bump ni publicación.<br>
**Último QA integrado:** `1db4f63b52af79247745b3a8a220fb728348218c`.

## Resultado

Se modernizó la superficie de PyQt6/QtWidgets manteniendo la lógica funcional y química existente: temas con tokens, recursos SVG/HiDPI, app bar y pestañas de documentos, rail con flyouts, panel lateral, barra de estado y CommandPalette (`Ctrl+P`). `Ctrl+K` conserva la acción Clean2D. La UI integrada recibió aprobación manual KDE/Wayland y no se inició una fase posterior.

### Decisiones que deben seguir vigentes

- No migrar a PySide/QML ni duplicar estado del canvas por razones visuales.
- Reutilizar QActions, docks/widgets, señales, shortcuts y handlers existentes; ocultar una superficie no significa eliminar su API.
- Mantener un único panel lateral, rail compacto y disposición compatible con 980×600; respetar HiDPI.
- Tratar Clean2D, ChemName, geometría de plantillas, persistencia `.cmsn` y química como campañas separadas.
- Conservar `docs/ui-modernization/pyqt6-spike/theme.py` como referencia de tokens porque la especificación global de tema aún la identifica normativamente. El ejecutable/prototipo visual se retiró; no es un entorno de desarrollo soportado.

## QA de cierre

- `compileall`: correcto; 1816 tests recolectados.
- Suite completa: 1760 passed, 55 skipped y un fallo conocido de CompChem async, sin regresión frente a baseline.
- Arquitectura: 269 passed; UI dirigida: 304 passed.
- Flatpak smoke instalado: PASS, incluidos assets UI y round-trip `.cmsn`.
- PyInstaller arrancó. El script del proyecto genera un ejecutable portable llamado `.AppImage`; no se afirmó formato AppImage Type 2 porque `appimagetool` no estaba instalado.
- Ruff scoped conserva un F401 histórico en un test Clean2D, fuera de alcance.

## Evidencia y memoria

- Capturas finales: [`after/README.md`](./after/README.md).
- Decisiones, desviaciones, intentos descartados, deudas y campañas Clean2D independientes: [`docs/history/CAMPAIGNS.md`](../history/CAMPAIGNS.md).
- Detalle de OpenSpec: [reporte QA de Fase 8](../../openspec/changes/ui-modernization-final-qa-release-2026-10-02/report.md) y cambios archivados bajo `openspec/changes/archive/2026-10-03-*`.
- La comparación homogénea before/after no existe: no se versionó captura PNG original de Fase 0, así que no se inventó un montaje.

Este documento es una ficha de estado, no una autorización para iniciar trabajo. Cualquier cambio nuevo requiere su propio alcance y OpenSpec.
