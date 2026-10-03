# Evidencia visual posterior — modernización UI

**Revisión:** 2 de octubre de 2026<br>
**Código:** `140a080c336515650bbeea0a4e6ead67a9999b23` + cierre documental de QA<br>
**Capturas:** PyQt6, plataforma `offscreen`, escala 100 %, ventana 1440×900.

| Captura | Contenido |
|---|---|
| `main-light-1440x900.png` | Ventana principal, tema claro. |
| `main-dark-1440x900.png` | Ventana principal, tema oscuro. |
| `toolrail-flyout.png` | Rail con flyout de enlaces abierto. |
| `sidepanel-templates.png` | Panel lateral en Plantillas. |
| `command-palette.png` | CommandPalette abierta (`Ctrl+P`). |

Las imágenes son capturas directas de la ventana y no se retocaron. El smoke también cubrió ambos temas a 980×600 y `QT_SCALE_FACTOR=2`: preview lógico de 88×56 puntos con backing store de 176×112 píxeles a DPR 2.

No hay captura PNG homogénea de Fase 0; no se fabricó montaje before/after. La aprobación manual KDE/Wayland de la UI integrada y capturas reales a DPR 2 se registran en [`visual-convergence/`](../visual-convergence/README.md). Para campañas, límites y resultado QA, ver [`docs/history/CAMPAIGNS.md`](../../history/CAMPAIGNS.md) y [reporte de Fase 8](../../../openspec/changes/ui-modernization-final-qa-release-2026-10-02/report.md).
