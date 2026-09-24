# Capturas Fase 3 — App bar y pestañas de documento

Capturas de referencia generadas **offscreen** desde la ventana real
(`ChemusonWindow`, `QT_QPA_PLATFORM=offscreen`), tras la Fase 3 de la
modernización de la UI (OpenSpec
`2026-09-24-modernize-ui-app-bar-document-tabs`).

## Escenarios

| Archivo | Escenario |
|---|---|
| `appbar-phase-full-light.png` | Ventana completa 1440×900, tema claro: 3 pestañas (una nombrada `Proyecto A.cmsn`, dos modificadas → punto de suciedad + sufijo `" *"` real), historial undo activo (botón deshacer habilitado), QMenuBar visible, toolbar clásico oculta, riel lateral y estado intactos. |
| `appbar-phase-full-dark.png` | Idem, tema oscuro: tinte de iconos y acento del tema oscuro, pestaña activa con subrayado de acento. |
| `appbar-phase-narrow-full-light.png` | Ventana 980×600 (comportamiento estrecho): las pestañas se eliden (`Sin tít…`), la píldora de búsqueda mantiene su ancho fijo de 250 px (igual que el spike), y lo esencial (pestañas, `+`, deshacer/rehacer, tema) sigue visible; el canvas permanece intacto y la altura de la barra es fija (54 px). |
| `appbar-phase-appbar-light.png` | Close-up de la `AppBar` en claro. |
| `appbar-phase-appbar-dark.png` | Close-up de la `AppBar` en oscuro. |

## Comparación con el spike aprobado

Referencia visual: `docs/ui-modernization/pyqt6-spike/app.py` (app bar de
54 px: marca con frasco + versión, pestañas de documento con acento en la
activa, botón `+` punteado, píldora de búsqueda "Buscar o ejecutar… Ctrl
K", deshacer/rehacer, separador, tema y preferencias).

- **Coincidencias**: orden de elementos, altura (54 px), glifos SVG
  `currentColor` a 1.75/24×24, píldora de 250 px con kbd hint, punto de
  suciedad de 7 px en el color de acento, pestaña activa con fondo
  `surface3` + subrayado de acento.
- **Diferencias intencionales**: la implementación usa `QWidget`
  nativos + tokens de diseño + QSS (no HTML/CSS literal); las pestañas son
  un espejo de `QTabWidget` real (no una lista estática), así que el
  contenido depende del estado real del documento (títulos, suciedad desde
  `QUndoStack`, orden), y la píldora de búsqueda no registra ningún atajo
  (Ctrl+K sigue siendo "Clean 2D full" hasta la Fase 6).

## Regenerar

Desde la raíz del repositorio (con el venv del proyecto):

```bash
QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
  python docs/ui-modernization/app-bar-phase-shots/make_shots.py \
  docs/ui-modernization/app-bar-phase-shots
```

El script solo lee la UI (no modifica producción); crea 3 pestañas, una
ruta ficticia de título (`/ficticio/Proyecto A.cmsn`, sin escribir
archivos) y un comando undo de no-op para provocar suciedad real.
