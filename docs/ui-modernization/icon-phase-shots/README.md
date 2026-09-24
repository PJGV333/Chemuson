# Capturas de referencia — Fase 2 (sistema de iconos SVG)

Ventana **real** de ChemUSON (no spike) con el set de iconos SVG de la
Fase 2 (`openspec/changes/2026-09-24-modernize-ui-svg-icons`), generada
offscreen (`QT_QPA_PLATFORM=offscreen`, `ChemusonWindow` + `grab()`).

| Archivo | Contenido |
|---|---|
| `icon-phase-light.png` / `icon-phase-dark.png` | Ventana completa (1280×860) en tema claro/oscuro. |
| `maintoolbar-*.png` | Barra superior (documentos, undo/redo, rotaciones, flip, clean, SMI). |
| `toolbar-draw-*.png` | Rail de dibujo (puntero, enlace, ancla ondulada, anillo, etiqueta, esfera, clean). |
| `toolbar-symbols-*.png` | Rail de símbolos (texto, corchetes, flecha, TLC, carga, energía, orbital). |
| `toolbar-text-*.png` | Barra de texto (16 px). |

## Cómo leer la comparación

- **Antes (Fase 1)**: `../foundation-shots/foundation-light.png` y
  `foundation-dark.png` (misma ventana, iconos QPainter de `gui/icons.py`).
- **Referencia visual aprobada**: spike PyQt6
  (`docs/ui-modernization/pyqt6-spike/`, commit `59e977d`) — el set SVG
  reproduce su lenguaje (24×24, trazo 1.75, esquinas redondeadas,
  `currentColor`).

Cambios esperados respecto a la Fase 1 (intencionales):
1. **Tinte**: de near-black/near-white a los tokens `icon`
   (`#475569` claro / `#C3CEDF` oscuro) — lenguaje monocromo del spike.
2. **undo/redo**: iconos SVG propios (se eliminó `QIcon.fromTheme` +
   pixel-loop de tint).
3. **Botón de orbitales**: ahora usa la geometría aprobada del spike
   (dos lóbulos), en lugar de la escalera MO de la versión QPainter.
4. **Átomos/esfera/SMI**: colores CPK (dominio), independientes del tema;
   la esfera lleva gradiente radial.

## Regenerar

```bash
PYTHONPATH=src uv run --no-project --offline \
  --python /home/unison-pjgv/Documentos/GitHub/Chemuson/.venv/bin/python \
  --with PyQt6 --with numpy --with Pillow --with rdkit --with certifi --with PyYAML \
  -- python docs/ui-modernization/icon-phase-shots/make_shots.py docs/ui-modernization/icon-phase-shots
```

(El script crea `ChemusonWindow`, aplica `light`→`dark` vía `_apply_theme()`
y captura la ventana y cada toolbar con `grab()`.)
