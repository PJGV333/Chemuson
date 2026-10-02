# Design: Pulido final de la UI moderna (Fase 7)

## Decisiones de diseño

### D1. Contraste: subir `text3` light, no rediseñar la paleta

El auditoría (Ornith, verificado por Qwen en `theme/tokens.py:61`) confirma que
solo el **light** falla: `text3 = #94A3B8` sobre `#F1F5F9` ≈ 2.3:1. El dark
(`#A5B4C7` sobre `#0B1120` ≈ 8.9:1) es correcto.

- `text3` light → `#5E6E82` (verificado: 4.76:1 sobre `#F1F5F9`, 5.21:1 sobre
  `#FFFFFF` y 4.98:1 sobre `#F8FAFC`, AA con margen en los tres fondos de
  texto secundario). `#64748B` (4.34:1 sobre `#F1F5F9`) se descarta por quedar
  bajo el umbral 4.5:1.
- `text2` light (`#475569`) y todos los tokens dark quedan intactos.
- El canvas sigue siendo hoja blanca (consistente con export PNG); el overlay
  de estado vacío sigue usando `text2`/`text3` por tokens, sin cambios
  adicionales.

### D2. Disabled: color/fondo con tokens soportados, `QAction` como fuente de verdad

El estado habilitado **no se copia**: sigue viniendo de la `QAction`
(`canUndoChanged` → `action_undo.setEnabled`, ya verificado en
`main_window.py`). El pulido es exclusivamente presentacional en QSS.

**Corrección post-push (gate manual Fase 7)**: Qt Style Sheets solo soporta
la propiedad `opacity` para `QToolTip`; no funciona en `QToolButton`,
`QPushButton`, `QLineEdit`, `QCheckBox`, etc. Se eliminan todas las reglas
`opacity:` añadidas en Fase 7 a widgets normales; no se introduce
`QGraphicsOpacityEffect` ni lógica nueva. Los estados `:disabled` quedan
inequívocos con propiedades soportadas:

- `QToolButton:disabled` → `color: text3` + fondo/borde transparentes (el
  hover ya no lo hace parecer activo porque `:disabled` aparece después de
  `:hover` en la hoja y gana).
- `QToolButton#railBtn:disabled` → `color: text3` + fondo `surface2` (antes
  transparente; con color de fondo el atenuado es inequívoco en light y dark).
- `QFrame[cls="flyItem"]:disabled` → fondo `surface3` (más marcado que el
  reposo `surface2`) + `border` + regla específica
  `QFrame[cls="flyItem"]:disabled #flyLbl { color: text3 }` (la etiqueta
  actual hereda `#flyLbl` text1 y parecía habilitada).
- `QFrame[cls="paletteItem"]:disabled` → fondo `surface2` + `color: text3`.
- `QPushButton:disabled` → fondo `borderStrong` + `text3` (ya existente).
- `QPushButton[flat="true"]:disabled` → `text3` + borde `border`.
- `QLineEdit:disabled` / `QSpinBox, QDoubleSpinBox:disabled` → fondo
  `surface2` + `text3` + `border-color: border`.
- `QCheckBox:disabled` / `QRadioButton:disabled` → `color: text3` (el
  indicador ya tiene `::indicator:disabled` con `surface2`/`border`).
- `QToolButton:disabled` de `get_tool_palette_stylesheet` y
  `#palette_grid QToolButton:disabled` → `text3` + fondo `surface3` + borde.
- Tests (`test_ui_polish.py`): se descarta la comprobación de opacidad;
  ahora se verifica (a) que **ninguna** regla de widget use `opacity` y (b)
  que cada selector `:disabled` auditable introduzca al menos un cambio
  visual soportado respecto al estado normal.
- No se introduce `QPainter` ni lógica de enabled en widgets: solo QSS.

### D3. Tooltips: convención `Nombre (Shortcut)` desde `QAction.shortcut()`

- El rail (`tool_rail.py`) tiene `QAction` históricos por botón
  (`spec.icon_action` / acción de flyout). Si `action.shortcut().toString()`
  no es vacío → tooltip = `"{nombre} ({shortcut})"`; si no → se conserva el
  tooltip del spec (las letras V/B/L/R/C/T/N/G/E/O son atajos contextuales
  propios del `ToolShortcutDispatcher`, no `QAction.shortcut()`, y se
  mantienen tal cual).
- `document_tabs.py`: el tooltip de "Nuevo" pasa a usar
  `action_new.shortcut().toString()` (`StandardKey.New` = Ctrl+N) en vez de la
  cadena hardcodeada.
- `SearchPill` (`app_bar.py`) ya usa "Buscar o ejecutar… (Ctrl+P)": se
  referencia al `QAction` real de apertura en el mismo tooltip (única
  fuente).
- No se duplica infraestructura: se reutiliza el patrón de
  `command_palette._action_tooltip` (que ya usa `entry.action.shortcut()`).

### D4. HiDPI: DPR-aware en bitmaps; QSS en px lógicos se conserva

Qt6 escala px lógicos de QSS/fuentes por devicePixelRatio; el riesgo real es
**rasterizar bitmaps en píxeles físicos**:

- `template_preview_icon` (`template_browser_service.py`): pixmap creado a
  `88 × dpr × 56 × dpr` con `painter.setDevicePixelRatio(dpr)`; margen, ancho
  de pen y tamaño de font escalados por `dpr`. El contenido (átomos/enlaces
  del grafo) **no cambia**: solo escala de render.
- `IconProvider` ya es DPR-aware (verificado: `_render_svg` usa
  `size * dpr`); `draw_glyph_icon` delega en `icon_dynamic` → OK.
- Tests razonables (no pixel-perfect): offscreen con `QT_SCALE_FACTOR=2`,
  la ventana se construye, rail/app bar/side panel existen, iconos
  no-nulos a dpr 1 y 2, y el pixmap de thumbnail tiene tamaño físico
  `88×56 × dpr`.
- **Limitación conocida (no resuelta en Fase 7)**: `_device_pixel_ratio`
  usa el DPR físico de la pantalla primaria (`screen.devicePixelRatio()`),
  consistente con `IconProvider`. El factor de escala manual de Qt
  (`QT_SCALE_FACTOR` en un monitor físico 1×) **no** se refleja en
  `devicePixelRatio()`, así que en ese caso extremo el thumbnail queda a
  `88×56` físicos y se interpola al escalar. El caso HiDPI real (pantallas
  de alto DPI físico) sí queda nítido. Se documenta como limitación, no
  como defecto del cambio.

### D4.1. `iconSize` del árbol de plantillas (corrección post-push, gate manual Fase 7)

`template_preview_icon()` genera pixmaps lógicos 88×56 DPR-aware, pero
`PlantillasDock.tree` no fijaba `iconSize`, así que `QTreeWidget` los reducía
al tamaño de icono por defecto (~16 px) y la previsualización quedaba
diminuta. Ajuste exclusivamente visual: `PlantillasDock` fija
`tree.setIconSize(QSize(88, 56))` (constante `_TEMPLATE_THUMB_SIZE`, coincide
con el tamaño lógico del pixmap). La altura de fila **no** requiere ajuste:
`QTreeWidget` dimensiona la fila al icono (fila hoja = 58 px ≥ 56 px del
icono, verificado offscreen a 100 % y 200 %; sin recorte). 88 px de ancho
cabe en el SidePanel de 340 px con margen para el texto del ítem. No se toca
grafos, molblocks, átomos, enlaces ni geometría química. Evidencia: capturas
antes/después del SidePanel (100 % y 200 %) en `evidence/templates_before_*.png` /
`templates_after_*.png` y montajes `templates_montage_*.png`.

### D5. Onboarding: overlay de shell, 3 pasos, persistido

- Nuevo `gui/onboarding.py` → `OnboardingOverlay(QWidget)`: hijo de la
  ventana, `WA_TranslucentBackground` + máscara con "agujero" que resalta la
  zona (rail / lienzo / panel derecho), tarjeta con título + texto +
  `Anterior`/`Siguiente`/`Cerrar` y checkbox "No mostrar de nuevo".
- Pasos fijos (texto del brief):
  1. Tool Rail — "Elige aquí las herramientas de dibujo y anotación."
  2. Canvas — "Dibuja, selecciona y edita tus estructuras en el lienzo."
  3. SidePanel — "Inspector, validación, propiedades, plantillas y
     apariencia están aquí."
- Inyección en `shell/assembly.py` tras montar el shell: si
  `application_settings().value("ui/onboarding/completed", False)` es falso →
  mostrar overlay. **Semántica de persistencia (corrección post-push, gate
  manual Fase 7)**: la señal del overlay es `finished(bool)`; el ensamblado
  fija `ui/onboarding/completed=True` **solo** si llega `True`:
  completar los 3 pasos → `True`; cerrar con "No volver a mostrar" marcado →
  `True`; cerrar sin marcar → `False` (NO se persiste; el onboarding se
  ofrece de nuevo en el siguiente arranque). Antes de la corrección
  `overlay.no_more` existía pero no participaba en la persistencia.
- Sin sistema paralelo: reutiliza `platform.settings` (mismo `QSettings` que
  el resto de la GUI). No toca `gui/canvas/` ni la escena: el overlay es un
  hijo de la ventana principal y posiciona el "agujero" sobre geometrías
  públicas (`tool_rail`, `centralWidget`, `side_panel`).
- `QuickStartDialog` legado se conserva (menú Ayuda) — no se elimina
  (fuera de alcance).

### D6. Plantillas en CommandPalette: refresh defensivo pequeño

Encadena la actualización existente (no se toca `TemplateLibrary` ni química):

```
mutación → template_controller → refresh_template_views
         → templates_menu (fuente de verdad)
_open_command_palette → _refresh_command_registry_templates (ya existe)
```

Gap: si la paleta está abierta y la biblioteca cambia, el registro no se
reconstruye. Arreglo: `_refresh_template_views` (main_window) también invoca
`_refresh_command_registry_templates()`. Deduplicación por identidad ya
existe; coste O(acciones) insignificante.

### D7. Iconos de formato de texto: 5 SVG + glifos B/I/U conservados

Se añaden `i-align-left.svg`, `i-align-center.svg`, `i-align-justify.svg`,
`i-subscript.svg`, `i-superscript.svg` (rejilla 24 px, trazo ~1.75,
`currentColor` vía tint del provider). La barra de formato real usa
`AlignLeft`/`AlignHCenter`/`AlignJustify` (sin `AlignRight`), por eso el
cuarto glifo es "justificar". Los botones de alineación/sub/sup
usan el provider; `B/I/U` (glifos) se conservan: son convención estándar de
editores de texto y no son "provisionales" en sentido de placeholder.

### D8. Referencias stale: solo donde la paleta es el sujeto

Se corrige "Ctrl+K" → "Ctrl+P" en `command_palette.py:1`, `qss.py:970`,
`tokens.py:167`, `main_window.py:312/316-317`, `app_bar.py:131`. Se
**conservan** `structure_actions.py` y `assembly.py` (ahí `Ctrl+K` es
correctamente Clean2D) y `tool_rail.py:396` (referencia a atajos de Clean2D).
`PLAN.md`/`mockup-ui.html` se conservan: son documentos históricos de diseño
(pre-Fase 6), no código activo.

### D9. Manual: reescritura de secciones, sin renumerar

`docs/MANUAL_USUARIO.md`:

- §2 → mapa moderno (AppBar, DocumentTabs, ToolRail + flyouts, SidePanel,
  texto de formato, lienzo, barra de estado, tema claro/oscuro, onboarding).
- §3 → referencias a "panel de herramientas" (no "barras izquierda/derecha").
- §6.5 → los paneles viven como tabs del SidePanel (el menú Ver sigue
  alternando su visibilidad, back-compat Fase 5).
- §10 → "Barra de aplicación" (tabs, pill `Ctrl+P`, undo/redo, tema) +
  subsección "Paleta de comandos (Ctrl+P)" con atajos reales
  (`Ctrl+K` Clean2D 1 paso, `Ctrl+Shift+K` publicación, `Ctrl+Alt+K`
  conformero).
- §11 → "Panel de herramientas (rail izquierdo + flyouts)" — mismo contenido
  de herramientas (sigue siendo exacto).
- §12 → "Herramientas de anotación y símbolos (flyouts del mismo rail)" —
  contenido conservado.
- §13 (barra de formato de texto) y §18 (paneles) se conservan (siguen
  exactos); §9 (Ayuda) gana una línea de la paleta de comandos como ruta
  rápida. No se renumeran §13–20.

## Riesgos y mitigación

- **QSS `opacity` sobre widgets**: Qt solo soporta `opacity` en `QToolTip`;
  no funciona en widgets normales. Corrección aplicada (gate manual Fase 7):
  las reglas `opacity:` añadidas en Fase 7 se eliminan y el estado
  `:disabled` se pinta solo con propiedades soportadas (color/fondo/borde de
  tokens); sin `QGraphicsOpacityEffect`.
- **Overlay sobre el canvas**: el overlay es hijo de la ventana (no de la
  escena); no se agregan items ni se toca hit-testing. El "agujero" solo
  posiciona una región transparente; el canvas sigue recibiendo eventos al
  cerrarse el paso.
- **Cambio de token `text3`**: afecta a todo texto terciario light; la
  verificación visual (screenshots) confirma que no se "ensucia" la UI.
- **DPR del thumbnail**: se prueban con `QT_SCALE_FACTOR=1` y `2`; sin
  pixel-perfect (solo tamaño físico y no-nulidad).
