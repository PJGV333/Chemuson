# Evidencia visual posterior — modernización UI

**Revisión:** 2 de octubre de 2026<br>
**Código:** `140a080c336515650bbeea0a4e6ead67a9999b23` + cierre documental en `release/ui-modernization-qa`<br>
**Captura:** PyQt6, plataforma `offscreen`, escala 100 %, ventana 1440×900.

| Captura | Contenido |
|---|---|
| `main-light-1440x900.png` | Ventana principal, tema claro. |
| `main-dark-1440x900.png` | Ventana principal, tema oscuro. |
| `toolrail-flyout.png` | Rail con el flyout de enlaces abierto. |
| `sidepanel-templates.png` | Panel lateral en Plantillas con plantillas incorporadas. |
| `command-palette.png` | Paleta de comandos abierta; `Ctrl+P`. |

Las cinco imágenes son capturas directas de la ventana, PNG de 1440×900; no se retocaron. El Smoke Qt también cubrió ambos temas en 980×600 y `QT_SCALE_FACTOR=2`. En HiDPI se comprobó un preview de 88×56 puntos lógicos con backing store de 176×112 píxeles a DPR 2.

La comparación before/after no incluye un montaje: el baseline de Fase 0 (`../baseline.md`) conserva descripción y datos del estado inicial, pero no hay capturas PNG originales versionadas que permitan una comparación visual homogénea. Estas imágenes automatizadas tampoco sustituyen la aprobación manual en KDE/Wayland, ya recibida para la interfaz integrada.
