# UI Tool Rail & Flyouts Specification

## Purpose

Define el rail de herramientas unificado y los flyouts de Chemuson (Fase 4):
una superficie visual vertical (~58 px) que reemplaza visualmente las dos
barras históricas de toolbars (`ChemusonToolbar`, `SymbolPaletteToolbar`),
reutilizando 1:1 sus `QAction`, `QActionGroup`, tool_ids, señales y
handlers, sin duplicar lógica de herramienta, sin tocar el contrato canvas
(`tool_id`/`set_current_tool`), sin tocar química/persistencia y sin
rediseñar los docks.

## ADDED Requirements

### Requirement: Rail de herramientas unificado (superficie visual, sin segunda lógica)

La aplicación SHALL exponer un `ToolRail` (`src/chemuson/gui/tool_rail.py`)
vertical de ancho fijo ~58 px en el área de toolbar izquierda, visible en la
ventana real, construido a partir de los toolbars existentes
(`ChemusonToolbar`, `SymbolPaletteToolbar`) como fuentes de verdad. El rail
SHALL delegar cada interacción a los `QAction`/callbacks de los toolbars
originales y SHALL NOT crear `QAction` nuevos de herramienta ni inventar
`tool_id` nuevos. Las toolbars históricas SHALL conservarse (no eliminarse)
y oculta, con su `QActionGroup` exclusivo intacto.

#### Scenario: La ventana real monta el rail y oculta las barras históricas
- **GIVEN** un `ChemusonWindow` instanciado (offscreen)
- **WHEN** se inspecciona la jerarquía de widgets
- **THEN** existe `window.tool_rail` visible de ancho fijo 58
- **AND** `window.toolbar` y `window.symbols_toolbar` existen pero no son visibles
- **AND** el `QMenuBar` es visible con los menús históricos
- **AND** ninguna `QMenu` de paleta histórica es visible.

#### Scenario: El rail expone todos los grupos de herramientas existentes
- **GIVEN** el `ChemusonWindow` montado
- **WHEN** se inspeccionan los botones del rail
- **THEN** existen botones para: seleccionar (V), lazo (A), enlace (B), cadena
  (L), anillo (R), átomo (C), coordinación, rotación 3D, texto (T), flecha
  (N), corchetes (G), símbolos, placas, diagramas de energía (E), orbitales
  (O), limpiar 2D, validar y numerar
- **AND** cada botón de paleta tiene un desplegable que abre su flyout
  **AND** cada `tool_id` del inventario del baseline (selección 2, enlace 11,
  anillo 11, átomo 10, corchetes 10, flecha 16, placas 2, símbolos 12,
  energía 8, orbitales 23) es alcanzable desde el rail.

### Requirement: Flyouts reutilizables que reemplazan visualmente los QMenu de paleta

El componente `Flyout` (`src/chemuson/gui/flyout.py`) SHALL ser un popover
frameless de ancho fijo 244 px (cabecera con título + `Kbd` "Esc",
cuadrícula de celdas icono 22 px + etiqueta con word-wrap 1–3 líneas, pie
opcional con hasta 3 botones de texto) que se muestre junto al botón del
rail que lo abre y se cierre con Esc, clic fuera o selección de celda. Las
celdas de los flyouts de paleta SHALL ejecutarse 1:1 con los callbacks de
los `QMenu` originales de los toolbars.

#### Scenario: Flyout de enlace migra la paleta de enlaces
- **GIVEN** la ventana montada
- **WHEN** se hace clic en el desplegable del botón de enlace (B)
- **THEN** se muestra un flyout con 11 celdas (mismas que el `QMenu`
  original del toolbar: sencillo, bold, doble, triple, aromático, wedge,
  hashed, wavy, flexible, interacción, coordinativo)
- **WHEN** se elige la celda "Wedge"
- **THEN** el toolbar original emite `bond_palette_changed` con
  `style==BondStyle.WEDGE` y `tool_changed` con `"tool_bond"`
- **AND** el flyout se cierra y el botón de enlace del rail queda activo.

#### Scenario: Flyout de anillo conserva la entrada de tamaño personalizado
- **GIVEN** la ventana montada
- **WHEN** se abre el flyout de anillo (R) y se pulsa su pie
  "Tamaño personalizado…"
- **THEN** se abre el `QInputDialog` del toolbar original (mismo callback).

#### Scenario: Flyout de átomo conserva la tabla periódica
- **GIVEN** la ventana montada
- **WHEN** se abre el flyout de átomo (C) y se pulsa su pie "Tabla
  periódica…"
- **THEN** el toolbar original emite `periodic_table_requested`.

#### Scenario: Flyout de energía migra presets y diálogos
- **GIVEN** la ventana montada
- **WHEN** se abre el flyout de diagramas de energía (E)
- **THEN** se muestran las 8 celdas de presets (mismos `tool_id`
  `tool_energy_diagram_*`) y un pie con 3 botones (Atómico…, MO…, Ligando…)
- **WHEN** se pulsa "MO…"
- **THEN** el toolbar original emite `diatomic_mo_diagram_requested`
- **AND** los 25 presets de los submenús originales siguen disponibles vía
  un menú que emite `electronic_diagram_preset_requested` con el nombre del
  preset.

#### Scenario: Flyout de orbitales reutiliza los iconos existentes
- **GIVEN** la ventana montada
- **WHEN** se abre el flyout de orbitales (O)
- **THEN** se muestran las 23 celdas de `ORBITAL_PALETTE_MODEL` con los
  `QIcon` de `draw_orbital_icon(kind)` (QPainter, sin cambio de aspecto)
- **AND** elegir una celda emite `tool_changed` con el `tool_id`
  `tool_orbital_{kind}` original.

#### Scenario: El flyout se cierra con Esc y clic fuera
- **GIVEN** un flyout abierto
- **WHEN** se pulsa Escape
- **THEN** el flyout se oculta y emite `closed`
- **WHEN** (en otra ejecución) se hace clic fuera del flyout
- **THEN** el flyout se oculta y emite `closed`.

### Requirement: Estado del rail sincronizado con el estado real de herramienta

El highlight activo del rail SHALL ser derivado (no la fuente de verdad) de
las señales `tool_changed` de los toolbars y del estado del canvas. La
fuente de verdad del estado de herramienta SHALL seguir siendo el canvas
(`state.active_tool`) + los toolbars (specs de paleta).

#### Scenario: Seleccionar herramienta actualiza el highlight del rail
- **GIVEN** la ventana montada con herramienta "tool_select"
- **WHEN** se selecciona "Wedge" en el flyout de enlace
- **THEN** el estado del canvas pasa a `tool_bond`
- **AND** el botón de enlace del rail queda marcado activo y el de selección
  deja de estarlo.

#### Scenario: Cambiar de pestaña limpia el highlight
- **GIVEN** la ventana con dos pestañas y una herramienta activa
- **WHEN** se cambia a la otra pestaña
- **THEN** `_clear_active_tool_selection` se ejecuta (canvas `tool_none`)
- **AND** el rail queda sin botones activos.

#### Scenario: No hay doble emisión de tool_changed
- **GIVEN** la ventana montada con un contador conectado a
  `toolbar.tool_changed`
- **WHEN** se selecciona una herramienta desde el rail
- **THEN** el contador recibe exactamente una emisión por selección.

### Requirement: Atajos contextuales de letra simple

La aplicación SHALL habilitar los atajos de letra simple `V, A, L, B, R, C, T,
N, G, E, O` (mapeo del diseño D4) mediante un dispatcher de event filter
contextual, SHALL NOT usar `QShortcut` para ellos, y SHALL not activarse con
modificadores, con diálogo modal activo o con el foco en un widget de
entrada de texto. `Ctrl+K` SHALL seguir asociado a `action_clean_2d_full`
(sin badge visible).

#### Scenario: Atajo selecciona herramienta con canvas activo
- **GIVEN** la ventana montada con el canvas con foco
- **WHEN** se pulsa `B` (sin modificadores)
- **THEN** el canvas cambia a `tool_bond` y el botón de enlace queda activo.

#### Scenario: Atajo inactivo sobre widget de entrada de texto
- **GIVEN** un `QLineEdit` (p. ej. de un diálogo abierto no modal) con foco
- **WHEN** se pulsa `R`
- **THEN** no se cambia la herramienta (la tecla se escribe).

#### Scenario: Atajo inactivo con modificadores o diálogo modal
- **GIVEN** la ventana montada
- **WHEN** se pulsa `Ctrl+B`
- **THEN** no se cambia la herramienta
- **AND** con un diálogo modal activo, `R` no cambia la herramienta.

#### Scenario: Ctrl+K sigue siendo Clean 2D full
- **GIVEN** la ventana montada
- **WHEN** se pulsa `Ctrl+K`
- **THEN** se dispara `action_clean_2d_full.trigger` (sin badge "Ctrl K"
  visible en la app bar).

### Requirement: Integración de tema y refresh de iconos

El rail y los flyouts SHALL usar el QSS de tokens (light/dark) y los
`METRICS` nuevos (`railW: 58`, `flyoutW: 244`), y el `refresh_icons()` del
rail SHALL re-leer los iconos de los `QAction` de los toolbars y
reconstruir las celdas de los flyouts tras cada cambio de tema.

#### Scenario: Tema light → dark → light en el rail
- **GIVEN** la ventana montada (light)
- **WHEN** se cambia a dark y luego a light
- **THEN** el rail y el flyout actualizan su estilo (palette/QSS de tokens)
- **AND** los iconos de los botones y celdas se regeneran sin pérdida de la
  selección actual (bond spec, elemento, etc. intactos).

#### Scenario: Métricas del rail
- **GIVEN** `chemuson.gui.theme.tokens.METRICS`
- **WHEN** se consultan
- **THEN** contienen `railW == 58` y `flyoutW == 244`.

## MODIFIED Requirements

(ninguna: los requisitos de las Fases 1–3 se conservan intactos; la
`main_toolbar` sigue oculta, la app bar y las pestañas no cambian).
