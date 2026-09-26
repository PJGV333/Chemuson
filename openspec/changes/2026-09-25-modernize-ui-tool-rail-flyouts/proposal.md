# Proposal: Rail de herramientas unificado + flyouts (Fase 4 del PLAN de modernización)

## Why

Tras la Fase 3 (app bar + pestañas de documento, polish aplicado), la ventana
sigue teniendo **dos barras verticales históricas** (`ChemusonToolbar`
izquierda, `SymbolPaletteToolbar` derecha) con botones-`QMenu` (cuadrícula de
`QToolButton` dentro de `QMenu`). El mockup aprobado
(`docs/ui-modernization/mockup-ui.html`) y el spike PyQt6 aprobado
(`docs/ui-modernization/pyqt6-spike/`, commit `59e977d`) definen el objetivo
de la Fase 4: **un único rail vertical (~58 px) a la izquierda** con botones
compactos (icono + kbd-hint) y **flyouts** (popover 244 px, cabecera,
cuadrícula de celdas icono+etiqueta, pie opcional) en lugar de los `QMenu`
históricos. Esta fase migra TODAS las paletas/herramientas existentes al rail
**sin perder funcionalidad** y **sin crear un segundo sistema de lógica de
herramientas**: el rail es una superficie visual que delega 1:1 a las
`QAction`, `QActionGroup`, tool_ids, señales y handlers existentes
(`toolbar.py`, canvas, main window).

## What Changes

- **`src/chemuson/gui/flyout.py`** (nuevo): `Flyout` — componente reutilizable
  (QFrame frameless, ancho fijo 244 px): cabecera (título + `Kbd` "Esc"),
  cuadrícula de `FlyoutCell` (icono 22 px + etiqueta con word-wrap 1–3
  líneas, `checked`/`active`), pie opcional (separador + etiqueta + hasta 3
  botones de texto), `show_near(anchor)`, cierre con Esc / clic fuera /
  selección, señal `closed`, callback `on_select`. Lenguaje visual del spike
  (padding 11 px, gap 6 px, radio 9 px, sombra `shadow2` de tokens).
- **`src/chemuson/gui/tool_rail.py`** (nuevo): `ToolRail` (QWidget vertical,
  ~58 px) + `ToolRailButton` (QToolButton autoRaise, icono 21 px + kbd-hint
  opcional + indicador de estado activo). Creado con los **toolbars
  existentes** como fuentes de verdad:
  - Botones de **paleta** (select, enlace, anillo, átomo, corchetes, flecha,
    texto, placas, símbolos, diagramas de energía, orbitales): un clic
    dispara el `QAction` del toolbar original (emite la herramienta actual);
    el desplegable (chevron / segundo clic) abre el `Flyout` construido
    **leyendo el `QMenu` original** del toolbar (mismas celdas, mismos
    callbacks vía `button.click()` de los `QToolButton` del menú) → cero
    duplicación de lógica; footer del flyout reutiliza los QActions de
    "Tamaño personalizado...", "Tabla periódica..." y los 3 diálogos de
    diagramas electrónicos (los 25 presets de los submenús se reabren en un
    `QMenu` conectado a la misma señal `electronic_diagram_preset_requested`).
  - Botones de **acción** (cadena, centro de coordinación, rotación 3D,
  limpiado 2D, validación, numeración): disparan el `QAction`/handler
  existentes.
  - Estado activo (highlight) derivado de las señales `tool_changed` de los
  toolbars + `clear_active()` desde `_clear_active_tool_selection` (cambio de
  pestaña); iconos de paleta re-lectos de los `QAction` del toolbar en
  `refresh_icons()`.
- **`src/chemuson/gui/shell/assembly.py`**: crear `ToolRail` tras los
  toolbars; añadirlo a `LeftToolBarArea` en un `QToolBar` mínimo wrapper
  (no movable/floatable); **ocultar** `self.toolbar` y
  `self.symbols_toolbar` (no eliminar: siguen poseyendo las `QAction`, el
  `QActionGroup` exclusivo, las señales y los callbacks); conectar
  `tool_changed` → `tool_rail.set_active_tool`; `_apply_theme` refresca el
  rail (tras los `refresh_icons` de los toolbars); atajos contextuales
  `V/A/L/B/R/C/T/N/G/E/O` vía `ToolShortcutDispatcher` (event filter de la
  ventana: letra simple, sin modificadores, foco no en
  `QLineEdit`/`QTextEdit`/`QPlainTextEdit`/`QComboBox`, sin diálogo modal
  activo) que disparan los mismos `QAction`/callbacks (sin `QShortcut` para
  evitar doble conexión).
- **`src/chemuson/gui/theme/tokens.py`**: métricas `railW: 58`, `flyoutW: 244`
  en `METRICS`.
- **`src/chemuson/gui/theme/qss.py`**: sección de QSS de tokens para
  `#toolRail`, `#railBtn` (+ `[railKind]`), `#railKbd`, `#railSep`, `#flyout`,
  `#flyLbl`, `#flyoutTitle`, `#flyFoot`, `#kbdPill`.
- **Tests** (`tests/test_ui_tool_rail.py`): inventario 1:1 (recuento de
  celdas por paleta vs menús originales), delegación (clic en celda del
  flyout → señales del toolbar original: `bond_palette_changed` con spec,
  `tool_changed`), footer periódico/size/ring, estado activo (tool_changed,
  clear en cambio de pestaña, exclusión del `QActionGroup`), atajos
  contextuales (fuego con canvas activo; inactivos en `QLineEdit` con foco,
  con modificadores, con diálogo modal), tema light→dark→light, toolbars
  antiguos ocultos y presentes, `QMenuBar` visible, sin `QMenu` de paleta
  visible.
- **Capturas** `docs/ui-modernization/tool-rail-phase-shots/` (light/dark,
  rail completo, flyouts de enlace/símbolos/orbitales abiertos) + script
  reproducible.
- **`architecture/modules.yml`**: registrar `tool_rail.py`, `flyout.py` en M08
  (paths + internal_api).

## No cambios

- **`toolbar.py`**: sin modificaciones de comportamiento (el rail solo lee sus
  `QMenu`/`QAction`/`QActionGroup`); `refresh_icons()` existente intacto.
- **Chemio / Clean2D / ChemName / persistencia**: intocados.
- **Docks** (right-side docks existentes): intactos (el mockup los muestra
  como side-panel; su migración es Fase 5, no se toca).
- **`main_toolbar`** (ya oculta en Fase 3): no se toca.
- **`orbitals.py`**: `draw_orbital_icon` (QPainter) se reutiliza tal cual en
  el flyout (residuo documentado; ver design.md §Orbitals).
- **`QActionGroup` exclusivo** y los `tool_ids`: sin ids nuevos de
  herramienta (el rail no inventa ids; solo espeja los existentes).
- **Fase 5**: no se inicia.

## Compatibilidad / rollback

- Rollback: revertir el commit; `toolbar.py`/canvas/ventana no cambian de
  comportamiento, así que al revertir el rail las barras históricas vuelven a
  mostrarse (el assembly oculta las toolbars tras montar el rail; en el
  rollback el assembly vuelve a mostrarlas).
- Las toolbars antiguas siguen siendo la fuente de verdad del estado
  (`_current_bond_spec`, `_current_element`, `current_*()`); el rail no
  guarda estado propio de herramienta (solo highlight derivado).
