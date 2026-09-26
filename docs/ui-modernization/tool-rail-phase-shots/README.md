# Capturas Fase 4 — Rail de herramientas unificado + flyouts

Capturas de referencia generadas **offscreen** desde la ventana real
(`ChemusonWindow`, `QT_QPA_PLATFORM=offscreen`), tras la Fase 4 de la
modernización de la UI (OpenSpec
`2026-09-25-modernize-ui-tool-rail-flyouts`).

## Escenarios

| Archivo | Escenario |
|---|---|
| `rail-phase-full-full-light.png` | Ventana completa 1440×900, tema claro: rail de 58 px en el borde izquierdo (18 botones: selección, lazo, enlace, cadena, anillo, átomo, coordinación, rotación 3D, texto, flechas, corchetes, símbolos, placas, energía, orbitales, limpiar 2D, validar, numerar), toolbars históricas ocultas, kbd hints (V/A/L/B/R/C/T/N/G/E/O), canvas intacto. |
| `rail-phase-full-full-dark.png` | Idem, tema oscuro: tokens oscuros del rail (`bg` #0F172A), acento de la herramienta activa. |
| `rail-phase-full-rail-light.png` | Close-up del rail en claro (58 px). |
| `rail-phase-full-rail-dark.png` | Close-up del rail en oscuro. |
| `rail-phase-bond-light.png` | Flyout de enlaces abierto (244 px): cabecera "ENLACES" + kbd `Esc`, 11 celdas (icono 22 px + etiqueta word-wrap) en 3 columnas; celda «Enlace doble» resaltada. |
| `rail-phase-bond-dark.png` | Idem, tema oscuro. |
| `rail-phase-energy-light.png` | Flyout de diagramas de energía abierto: 8 celdas de preset (2 columnas) + pie con los 2 submenús originales como botones de texto ("Diagrams ▾", "Presets ▾"). |
| `rail-phase-energy-dark.png` | Idem, tema oscuro. |

## Regeneración

    QT_QPA_PLATFORM=offscreen PYTHONPATH=src \
        .venv/bin/python docs/ui-modernization/tool-rail-phase-shots/make_shots.py \
        docs/ui-modernization/tool-rail-phase-shots

## Comparación con el mockup / spike

Referencia visual: `docs/ui-modernization/mockup-ui.html` (`.rail` 58 px,
`.flyout` 244 px con padding 11 px, gap 6 px, radio 9 px, celda con icono
22 px + etiqueta de 1–3 líneas) y el spike PyQt6 aprobado
(`docs/ui-modernization/pyqt6-spike/widgets.py`: `RailButton`,
`FlyoutCell`, `Kbd`, `SearchPill`).

- **Coincidencias**: ancho del rail (58 px) y del flyout (244 px), botones
  de 44 px con kbd hint en la esquina, cabecera del flyout con título en
  mayúsculas + kbd `Esc`, celdas con icono 22 px y etiqueta con word-wrap
  de 1–3 líneas, separadores de sección, pie con botones de texto
  (1–3), estados hover/activo con los tokens de acento.
- **Diferencias intencionales**: la implementación usa `QWidget` nativos +
  tokens de diseño + QSS (no HTML/CSS literal); el flyout se abre con
  clic derecho en el botón de paleta (el clic izquierdo activa la
  herramienta *actual*, igual que la barra histórica) y se cierra con
  `Esc`, clic fuera o selección de celda; la sombra del flyout se aplica
  con `QGraphicsDropShadowEffect` (el QSS de Qt no soporta `box-shadow`);
  el estado activo del rail es *derivado* de las señales `tool_changed`
  de los toolbars (el canvas + los toolbars siguen siendo la fuente de
  verdad; no hay `QAction` ni `tool_id` inventados).

## Notas de paridad 1:1 (baseline)

- 18 botones del rail → las mismas `QAction`/callbacks de los toolbars
  históricos (`ChemusonToolbar` + `SymbolPaletteToolbar`), que quedan
  **ocultos, no eliminados** (propietarios del `QActionGroup` exclusivo,
  del estado de paleta y de las señales).
- Celdas de flyout por paleta (equivalentes a los `QMenu` históricos):
  selección 2, enlaces 11, anillos 11, átomos 10, flechas 16, corchetes 10,
  símbolos 12, placas 2, energía 8, orbitales 23 (los 5 slots vacíos del
  grid de 28 del menú de orbitales no representan una función y se omiten).
- Pies de flyout: anillo → "Tamaño personalizado…" (`QAction` original),
  átomo → "Tabla periódica…" (`QAction` original), energía → 2 submenús
  originales re-abiertos vía `popup()` (mismas `QAction` de preset), texto
  → submenú "Color de etiquetas".
- Atajos de letra simple contextuales: V/A/L/B/R/C/T/N/G/E/O vía
  `ToolShortcutDispatcher` (event filter; **sin `QShortcut`**): se suprimen
  con modificadores, diálogo modal activo o foco en `QLineEdit`/
  `QTextEdit`/`QPlainTextEdit`/`QComboBox`/`QAbstractSpinBox`.
- Residuo documentado: los iconos de orbitales reutilizan el `QPainter`
  existente (`draw_orbital_icon` en `gui/orbitals.py`); la migración a SVG
  se pospone (cambio de comportamiento visual de 23 iconos, fuera del
  alcance quirúrgico de la fase).
