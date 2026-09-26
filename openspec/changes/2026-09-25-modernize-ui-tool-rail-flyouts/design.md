# Design: Rail de herramientas unificado + flyouts (Fase 4)

## Contexto y auditoría previa

Ver `baseline.md` (inventario 1:1 de tool_ids, acciones, señales, atajos y
decisión de `orbitals.py`). Resumen de los hechos que dirigen el diseño:

- `ChemusonToolbar` (izquierda) y `SymbolPaletteToolbar` (derecha) poseen el
  `QActionGroup` exclusivo compartido, las `QAction` de cada botón, los
  callbacks de selección de paleta (`_select_bond_palette`,
  `_select_ring_palette`, `_select_element_palette`, `_select_bracket_tool`,
  `_select_arrow_tool`, `_select_plate_tool`, `_select_symbol_tool`,
  `_select_energy_diagram_tool`, `_select_orbital_tool`,
  `_select_selection_palette`), los estados (`_current_bond_spec`,
  `_current_ring_spec`, `_current_element`, `current_*()`), y emiten todas las
  señales (`tool_changed`, `bond/ring/element_palette_changed`,
  `periodic_table_requested`, 4 señales de diagramas electrónicos).
- El canvas consume los `tool_id` en `set_current_tool`
  (`canvas_selection_input.py`); la ventana actualiza estado y barra de
  estado (`_on_tool_changed`, `_handle_*_palette`, `_update_status`);
  `_clear_active_tool_selection` limpia al cambiar de pestaña.
- Ningún atajo de letra simple existe hoy (todos los atajos llevan
  modificadores o son F-keys) → `V/A/L/B/R/C/T/N/G/E/O` son libres.
- `draw_orbital_icon` usa `QPainter` (renderer propio con gradientes y
  `QPainterPath`); migrarlo a SVG cambiaría el aspecto → se reutiliza tal
  cual (residuo documentado).

## Decisiones de diseño

### D1 — El rail es superficie visual; la lógica vive en los toolbars

`ToolRail(toolbar: ChemusonToolbar, symbols_toolbar: SymbolPaletteToolbar)`
no crea `QAction` nuevos de herramienta ni duplica specs:

- **Botón de acción simple** (cadena, coordinación, 3D, limpiado, validar,
  numerar): el clic dispara la `QAction`/handler existente
  (`chain_action.trigger()`, `coord_action.trigger()`,
  `rotate_3d_precise_action.trigger()`, `window.action_clean_2d_full.trigger()`,
  `action_validate_structure.trigger()`, `action_numbering_recalculate.trigger()`).
- **Botón de paleta**: un clic dispara el `QAction` original del toolbar
  (`select_action`, `bond_action`, `ring_action`, `label_action`,
  `bracket_action`, `annotation_action`, `plate_action`, `symbol_action`,
  `energy_diagram_action`, `orbital_action`, `text_action`) → el toolbar emite
  su señal habitual con la selección *actual*. El desplegable abre el
  `Flyout` construido **leyendo el `QMenu` original** del botón
  (`button.menu()`): se recorre la `QWidgetAction` (contenedor con
  `QGridLayout`) y se toma cada `QToolButton` (icono, tooltip, su callback
  conectado) + las `QAction` no-widget del menú (footers: "Tamaño
  personalizado...", "Tabla periódica...", submenús de diagramas
  electrónicos). El clic en una celda del flyout llama `button.click()` del
  `QToolButton` del menú original → el callback del toolbar original se
  ejecuta (actualiza iconos/tooltips/estado del toolbar y emite las señales
  originales). **Cero lógica de herramienta nueva; cero `tool_id` nuevo.**
- El `QActionGroup` exclusivo se conserva intacto en el toolbar original: el
  estado checkado real de la herramienta vive ahí; el highlight del rail es
  derivado (D3), no es la fuente de verdad.

### D2 — `Flyout` (componente reutilizable, `gui/flyout.py`)

`QFrame` frameless, `setFixedWidth(METRICS["flyoutW"])` (244 px), objeto
hijo de la ventana (no modal; cierra con Esc, clic fuera o selección —
event filter de la ventana; `closeOnSelect=True` por defecto).

- `populate(title: str, items: list[FlyoutItem], columns: int,
  active_id: str | None, on_select: Callable[[str], None],
  footer: FlyoutFooter | None) -> None`
  - `FlyoutItem(id, icon: QIcon, label: str, enabled: bool = True,
    tooltip: str = "")`
  - `FlyoutFooter(text: str, buttons: list[tuple[str, Callable[[], None]]])`
    (1–3 botones de texto)
- `FlyoutCell` (equivalente del spike): icono 22 px + `QLabel` word-wrap
  (10 px, 1–3 líneas, wrap manual con `QFontMetrics`), `minimumHeight=56`,
  `active`/`checked` vía property + `unpolish/polish`, hover/active con
  `accentSoft`/`accentBorder`.
- `show_near(anchor: QWidget, parent: QWidget)`: posición `x =
  anchor.mapTo(parent, (width, height//2))` + 8 px, centrado verticalmente
  sobre el botón, clamp al rect de la ventana (si no cabe a la derecha,
  aparece a la izquierda).
- Señal `closed()`; callback `on_select(item_id)` (al elegir una celda:
  invoca + oculta).
- Pie: separador (1 px `border`), etiqueta (`text3`) y botones de texto
  (`#flyFoot`); los botones no cierran el flyout salvo que el callback lo
  indique (`closeAfter: bool` por botón, default True).

### D3 — Sincronización estado rail ↔ estado real de herramienta

- El rail se conecta (en el assembly) a `toolbar.tool_changed` y
  `symbols_toolbar.tool_changed` → `tool_rail.set_active_tool(tool_id)`.
  Mismo mecanismo que `_update_status` (no se añade señal nueva ni estado
  nuevo).
- `set_active_tool` resuelve el grupo de rail para cada id:
  - `tool_select`/`tool_select_lasso` → botón de selección (el icono se
    re-lee de `toolbar.select_action.icon()` para distinguir pointer/lasso).
  - `tool_bond` → botón de enlace (icono de `bond_action`).
  - `tool_ring` → anillo; `tool_atom` → átomo (icono de `label_action`,
    que ya refleja el elemento activo vía `_set_element_palette_ui`).
  - `tool_brackets_*` → botón de corchetes; `tool_arrow_*` → flecha;
    `tool_tlc`/`tool_electrophoresis` → placas; `tool_charge*`/
    `tool_symbol_*` → símbolos; `tool_energy_diagram_*` → energía;
    `tool_orbital_*` → orbitales; `tool_text` → texto; `tool_chain`,
    `tool_coordination_center`, `tool_rotate_3d_precise` → flash
    transitorio (no es un "tool" persistente: se marca activo durante 400 ms
    y se desmarca; el estado real lo gestiona el canvas).
  - `tool_none`/otro → sin highlight.
- Iconos de paleta: `refresh_icons()` del rail re-lee `icon()` de los
  `QAction` de los toolbars y reconstruye las celdas de los flyouts (llamado
  en `_apply_theme` **después** de los `refresh_icons()` de los toolbars,
  que son los que regeneran los iconos).
- Cambio de pestaña: `_clear_active_tool_selection` llama
  `tool_rail.clear_active()` (una línea; el canvas ya recibe `tool_none`).
- No hay doble emisión: el rail nunca emite `tool_changed`; solo dispara
  `QAction.trigger()`/callbacks de los toolbars, que son la única fuente de
  señales.

### D4 — Atajos contextuales (sin `QShortcut`)

`ToolShortcutDispatcher` (QObject, event filter de la ventana, creado en el
assembly):

- Mapa (auditado en baseline.md): `V→select`, `A→lasso`, `L→cadena`,
  `B→enlace`, `R→anillo`, `C→átomo`, `T→texto`, `N→flecha`, `G→corchetes`,
  `E→energía`, `O→orbitales`. Cada entrada dispara el **mismo**
  `QAction`/callback que el clic en el rail (D1).
- Filtro (KeyPress, solo si la ventana es `activeWindow`):
  1. `event.modifiers() == 0` (sin Ctrl/Alt/Shift/Meta).
  2. `QApplication.activeModalWidget() is None`.
  3. El widget con foco (o ancestro en foco) no es `QLineEdit`, `QTextEdit`,
     `QPlainTextEdit`, `QComboBox` ni `QAbstractSpinBox` (no robar la
     escritura).
  4. La tecla está en el mapa → dispara el callback y `event.accept()`.
- No se usan `QShortcut` (evita conflictos de contexto y doble conexión); no
  hay atajo para `Ctrl+K` (sigue en `action_clean_2d_full`, sin badge,
  Fase 6 según polish Fase 3).

### D5 — Layout y montaje (assembly)

```
QMainWindow
├── QMenuBar (visible, histórico)
├── toolBar(LeftToolBarArea):  QToolBar wrapper (no movable) → ToolRail (~58 px)
│                              (self.toolbar y self.symbols_toolbar: setVisible(False))
├── toolBar(TopToolBarArea):   main_toolbar (oculta, Fase 3)
└── central QWidget (VBox)
    ├── AppBar (54 px)
    ├── TextFormatToolbar
    └── QTabWidget (tabs)
```

- `ToolRail` se añade en un `QToolBar` wrapper mínimo
  (`setMovable(False)`, `setFloatable(False)`, sin título; `addWidget(rail)`);
  el wrapper se estiliza con el QSS del rail (mismo `surface`, sin chrome).
- Orden de montaje: toolbars → rail (necesita las `QMenu`/`QAction` ya
  construidas) → ocultar toolbars → conexiones de señal.
- `TextFormatToolbar`, docks, `main_toolbar`: no se tocan.
- `orbitals.py`: `draw_orbital_icon(kind)` (QPainter) se reutiliza en las
  celdas del flyout de orbitales y en el botón del rail (icono de
  `orbital_action`). **Residuo documentado**: la migración a la
  infraestructura SVG de la Fase 2 queda pendiente (cambiar el renderer
  implicaría redibujar 23 orbitales y riesgo de drift visual; fuera de esta
  fase).

### D6 — Orden y contenido del rail (mockup + parity 1:1)

| # | Rail button (kbd) | Owner (toolbar) | Flyout |
|---|---|---|---|
| 1 | Seleccionar (V) | `toolbar.select_action` (emite la selección actual) | 2 celdas: pointer/lasso (leído del menú original) |
| 2 | Lazo (A) | callback `toolbar._select_selection_palette("tool_select_lasso", …)` | — |
| — | separador | | |
| 3 | Enlace (B) | `toolbar.bond_action` | 11 celdas (leído del menú; callback `button.click()`) |
| 4 | Cadena (L) | `toolbar.chain_action.trigger()` | — |
| 5 | Anillo (R) | `toolbar.ring_action` | 11 celdas + footer "Tamaño personalizado…" (QAction del menú) |
| 6 | Átomo (C) | `toolbar.label_action` | 10 celdas + footer "Tabla periódica…" (emite `periodic_table_requested`) |
| 7 | Coordinación | `toolbar.coord_action.trigger()` | — |
| 8 | Rotación 3D | `toolbar.rotate_3d_precise_action.trigger()` | — |
| — | separador | | |
| 9 | Texto (T) | `symbols_toolbar.text_action` | celdas = QActions del menú de formato (11) + colores (2) |
| 10 | Flecha (N) | `symbols_toolbar.annotation_action` | 16 celdas |
| 11 | Corchetes (G) | `symbols_toolbar.bracket_action` | 10 celdas |
| — | separador | | |
| 12 | Símbolos | `symbols_toolbar.symbol_action` | 12 celdas |
| 13 | Placas | `symbols_toolbar.plate_action` | 2 celdas (TLC, Gel) |
| — | separador | | |
| 14 | Energía (E) | `symbols_toolbar.energy_diagram_action` | 8 celdas + footer 3 botones (Atómico…, MO…, Ligando… → QActions de diálogos) + menú de 25 presets (señal `electronic_diagram_preset_requested`) |
| 15 | Orbitales (O) | `symbols_toolbar.orbital_action` | 23 celdas (grid 4×7, como `ORBITAL_PALETTE_MODEL`) |
| — | separador | | |
| 16 | Limpiar 2D | `window.action_clean_2d_full.trigger()` | — |
| 17 | Validar | `window.action_validate_structure.trigger()` | — |
| 18 | Numerar | `window.action_numbering_recalculate.trigger()` | — |

Contadores verificados por test: 2 / 11 / 11(+1) / 10(+1) / 16 / 10 / 12 / 2 /
8(+3) / 23 (+ 13 QActions de formato/texto del menú de texto).

### D7 — QSS y tokens

- `METRICS`: `railW: 58`, `flyoutW: 244`.
- QSS (tokens): `#toolRail { background: surface; border-right: 1px solid
  border; }`, `#railBtn` (autoRaise, padding 8, radio `radiusBtn`, hover
  `surface2`, active `accentSoft` + border `accentBorder` + color
  `accentStrong`), `#railKbd` (9 px, `text3`), `#railSep`, `#flyout`
  (surface, border, radio 10, sombra `shadow2` vía `QGraphicsDropShadowEffect`),
  `#flyLbl` (10 px, `text2`), `#flyoutTitle` (11 px bold, `text3`,
  uppercase), `#flyFoot` (11 px, `text2`, hover `surface2`).
- El wrapper `QToolBar` del rail hereda el QSS principal (tokens) — no se le
  aplica `TOOL_PALETTE_STYLESHEET` (ese es de las barras históricas).
- `refresh_icons()` (rail): re-lee iconos de los `QAction` y reconstruye
  celdas; `_apply_theme` lo invoca tras los de los toolbars.

## Riesgos y mitigaciones

- **Doble conexión de señales**: el rail no conecta señales nuevas de
  `tool_changed` a handlers de lógica; solo `set_active_tool` (visual). Test:
  un `tool_changed` actualiza una vez el estado del canvas.
- **Introspección de `QMenu`**: si un menú de paleta cambia de estructura,
  el rail lo detecta (assert de contadores en tests; en runtime, flyout
  vacío + warning log). Los menús actuales son estables (build una vez +
  `refresh_icons` los reconstruye con la misma estructura).
- **Atajos robando escritura**: filtro D4 + tests con `QLineEdit`/`QTextEdit`
  con foco, modificadores y diálogo modal.
- **Iconos QPainter de orbitales en flyout**: aspecto idéntico al toolbar
  actual (mismo `QIcon`); solo cambia el contenedor.
- **`QActionGroup` exclusivo**: intacto en el toolbar; el rail no lo toca
  (highlight derivado). Test: seleccionar paleta deja el check en la
  `QAction` del toolbar.
