# Design: Barra de aplicación y pestañas de documento (Fase 3)

## Contexto y auditoría previa

### Mapa actual: acción → QAction → superficie visual hoy

| Acción existente | QAction (dónde se crea) | Handler (main_window) | Superficie hoy | Superficie Fase 3 |
|---|---|---|---|---|
| Nuevo | `action_new` (Ctrl+N) — `actions/file_actions.py` | `_on_file_new` | menú Archivo + main_toolbar | **app bar `+`** (`action_new.trigger()`) + menú + Ctrl+N |
| Abrir | `action_open` (Ctrl+O) — `actions/file_actions.py` | `_on_file_open` | menú Archivo + main_toolbar | menú + Ctrl+O |
| Guardar | `action_save` (Ctrl+S, Ctrl+G) — `actions/file_actions.py` | `_on_file_save` | menú Archivo + main_toolbar | menú + Ctrl+S/Ctrl+G |
| Deshacer | `action_undo` (Ctrl+Z) — `actions/edit_actions.py` | `_on_undo` (estado real vía `canvas.undo_stack.canUndoChanged → action_undo.setEnabled`) | menú Editar + main_toolbar | **app bar undo** (mismo QAction) + menú + Ctrl+Z |
| Rehacer | `action_redo` (atajo estándar de Redo de la plataforma) — `actions/edit_actions.py` | `_on_redo` (id.) | menú Editar + main_toolbar | **app bar redo** (mismo QAction) + menú + atajo estándar de Redo de la plataforma |
| Preferencias | `action_preferences` — `main_window_ui_builder.build_menu_bar` | `_on_preferences` | menú Editar | **app bar settings** (mismo QAction) + menú |
| Tema | (no había QAction) — método `ChemusonWindow.toggle_theme` | `toggle_theme(checked)` | diálogo Preferencias (combo tema) | **app bar theme** (nueva `action_theme_toggle` checkable → handler existente) |
| Copiar/Cortar/Pegar/Duplicar | `action_copy`/`_cut`/`_paste`/`_duplicate` (Ctrl+C/X/V/D) | `_on_copy`… | menú Editar + main_toolbar (aux) | menú + atajos |
| Rotar/Flip | `action_rotate_left/right`, `action_flip_horizontal/vertical` | `_on_rotate…` | menú Editar→Rotar + main_toolbar | menú + atajos (Ctrl+Alt+Left/Right/I/A son de ramas; los de selección no tienen atajo) |
| Clean 2D | `action_clean_2d` (+ full Ctrl+K, publication Ctrl+Shift+K, propose Ctrl+Alt+K) | … | menú Estructura + main_toolbar | menú + atajos |
| SMILES | `action_draw_smiles` / `action_import_smiles` | `_on_import_smiles` | menú Estructura + main_toolbar | menú |

**Conclusión**: al ocultar `main_toolbar`, TODA acción sigue accesible por
menú, atajo y/o app bar (criterio de parada verificado por test).

### Estado real de documentos (fuente única de verdad)

- `QTabWidget self.tabs` (central) + `CanvasTabManager`
  (`tab_manager.py`): `create_document_tab` / `discard_canvas` /
  `close_canvas_tab` / `set_canvas_file_path` / `update_tab_title`.
- **Suciedad real**: `canvas.undo_stack.isClean()` (QUndoStack).
  `update_tab_title` pinta `" *"` en el texto de la pestaña nativa.
- Títulos base: `tab_manager.tab_titles[canvas]`; rutas:
  `tab_manager.file_paths[canvas]`.
- Señales útiles: `tabs.currentChanged`, `tabs.tabCloseRequested`,
  el callback `on_tab_updated` del `CanvasTabManager` (invocado al final de `update_tab_title`, tras el `setTabText` de `" *"`),
  `tabs.tabBar().tabMoved` (drag nativo). **No existe señal de
  add/remove** en `QTabWidget` → se añade un observer opcional en
  `CanvasTabManager` (D2).
- Creación de pestañas siempre pasa por `CanvasTabManager.create_document_tab`
  (assembly, `FileController.open_file_path`, `RecoveryController`,
  reemplazo tras cerrar); remoción siempre por `discard_canvas`
  (`_close_canvas_tab`, `FileController` en error de apertura,
  `RecoveryController`). Un observer en esas dos cubre **todos** los caminos.
- `Ctrl+K` está ocupado hoy por `action_clean_2d_full` → la píldora de
  búsqueda **no registra atajo** en esta fase (D6).

## Decisiones

### D1 — AppBar como `QFrame` en wrapper central, no `QToolBar`

`QMainWindow` solo acepta toolbars/menubar en su franja superior; un
`QToolBar` introduce comportamiento de desbordamiento (botón `>>`) y
manijas que no queremos en la barra de marca. Decisión: el widget central
pasa de ser directamente `self.tabs` a un wrapper `QWidget` con VBox
`[app_bar, self.tabs]` (margen 0, spacing 0, como el spike). El `QTabWidget`
sigue existiendo exactamente igual (mismo objeto, mismas señales, mismas
páginas); solo cambia su contenedor. Único uso de `centralWidget()` en el
repo es el propio assembly → sin impacto.

### D2 — `DocumentTabBar` = espejo de solo-lectura del `QTabWidget`

La nueva barra de pestañas **no guarda estado de documento**: títulos desde
`tab_manager.tab_titles`, suciedad desde `canvas.undo_stack.isClean()`,
órdenes e índices desde el `QTabWidget`. El `QTabWidget` oculta su tab bar
nativa (`.tabBar().hide()`); las páginas siguen controladas por
`currentIndex`. Espejo = `QTabBar` real (elide right, scroll buttons,
`movable`) + por pestaña un widget derecho (punto de suciedad 7 px `accent`
+ botón cerrar 18 px) + corner widget `+` (28 px). Señales del espejo:
`tabBarClicked → tabActivated(int)`, `tabCloseRequested(int)` (nativa),
`tabMoved(int,int)` (nativa), `newDocumentRequested()`.

Protocolo de sincronización (una sola dirección de datos hacia el espejo):

1. **Estructural (add/remove)**: `CanvasTabManager` recibe kwarg opcional
   `on_change: Optional[Callable[[CanvasTabManager], None]] = None` y lo
   invoca al final de `create_document_tab` y `discard_canvas`.
   Por defecto `None` → cero cambio de comportamiento para callers/tests
   existentes. La ventana pasa `self._on_document_tabs_changed` →
   resincronización **completa** del espejo (reconstruye pestañas: barata,
   idempotente, sin deriva).
2. **Texto/suciedad**: callback `on_tab_updated(canvas)` del `CanvasTabManager` (fue de `update_tab_title`) → slot que lee el
   estado real (`tab_titles.get(canvas)`, `not canvas.undo_stack.isClean()`)
   y actualiza esa pestaña del espejo. (Es el mismo estado que usa
   `update_tab_title`.)
3. **Selección**: el espejo emite `tabActivated(i)` →
   `tabs.setCurrentIndex(i)` → `currentChanged` → `_on_tab_changed`
   (handler existente: activa canvas, sync de acciones) →
   `app_bar.select(i)` (idempotente). Sin bucles (solo `tabBarClicked`
   dispara).
4. **Cerrar**: botón cerrar → `_on_tab_close_requested(i)` (handler
   existente: confirmación de descarte, reemplazo si queda 0, reactivación).
5. **Nuevo**: `+` → `action_new.trigger()` (QAction existente →
   `_on_file_new` → observer → resync).
6. **Reordenar**: drag en el espejo (`tabMoved`) →
   `tabs.tabBar().moveTab(from, to)` (el `QTabWidget` reordena sus páginas
   y emite `currentChanged` si aplica).

### D3 — `QMenuBar` visible y funcional; `main_toolbar` clásica oculta

Migración conservadora: el `QMenuBar` (Archivo/Editar/Ver/Estructura/
Reacción/Ayuda) **permanece visible y accesible** (Alt nativo, atajos
globales intactos). La franja superior clásica que la app bar reemplaza es la
`main_toolbar`: se oculta con `setVisible(False)` — el objeto, sus QAction y
`refresh_main_toolbar_icons` se conservan (compatibilidad; la Fase 4
reestructura el shell). Límite conocido documentado: el item de menú "Mostrar
copiar/pegar/zoom en barra superior" sigue operando sobre la toolbar oculta
(sin efecto visible hasta que se vuelva a mostrar).

### D4 — Reutilización de QAction (sin duplicar handlers ni atajos)

- Undo/redo/settings: `QToolButton.setDefaultAction(qaction)` → el botón
  hereda icono (ya puesto por `refresh_main_toolbar_icons`), tooltip
  (texto + atajo) y **estado enabled real** (Qt atenúa el icono al
  deshabilitar).
- `+`: `action_new.trigger()` (no se crea segundo "nuevo").
- Tema: **una** nueva QAction checkable `action_theme_toggle`
  ("Cambiar tema claro/oscuro") conectada al handler existente
  `toggle_theme` (sin duplicar lógica). Se crea en
  `create_local_actions`; su `checked` se sincroniza con
  `current_theme == "dark"` en `_apply_theme` (blockSignals). Icono:
  `moon` en claro, `sun` en oscuro (indicación del destino).
- `action_preferences` se mueve de `build_menu_bar` a
  `create_local_actions` (mismo objeto/conexión/menú) para que exista cuando
  se construye la app bar (orden de assembly).
- **Ningún atajo nuevo**: la píldora de búsqueda no registra shortcut
  (Ctrl+K pertenece a `action_clean_2d_full`).

### D5 — Píldora de búsqueda: placeholder de la Fase 6

`SearchPill` (QFrame): icono `search` (14 px, `text3`), texto "Buscar o
ejecutar…", kbd "Ctrl K" como pista visual, `clicked` emite señal
`searchActivated()` sin conectar a nada en esta fase (la command palette es
Fase 6 y traerá su atajo y su lógica). Tooltip: "(próximamente)".

### D6 — Iconos: 8 SVG nuevos del spike + `IconProvider` existente

`plus`, `search`, `moon`, `sun`, `sliders`, `flask`, `x`, `doc` (24×24,
`currentColor`, trazo 1.75; geometrías del spike aprobado). El inventario
pasa de 55 a 63 archivos. Tintes: marca `accent`; botones de app bar
`icon`; pill `text3`; tab icon `text2`/`text3`; close `text3`. Todo vía
`IconProvider` (caché por nombre+color+tamaño; HiDPI; sin pixel loops). La
marca usa `pixmap()` (QLabel), nunca `icon()` (qFatal QLabel/QIcon).

### D7 — QSS de tokens (Fase 1) extendida

Sección nueva en `get_main_stylesheet` (solo tokens, sin `box-shadow`/
`transitions`): `#app_bar` (surface + border-bottom, 54 px por
`METRICS["appbarH"]`), `#appBrandName` (14 px/700), `#appVersionPill`,
`#docTabs::tab` (8 px radio, selected = surface3 + borde inferior 2 px
`accent`, hover = surface2), `#tabNewBtn` (dashed, hover accent),
`[tabClose]` (hover border), `#searchPill`/`#searchPillTxt`/`#kbdK`,
`[appBarBtn]` (hover/pressed surface2/3), `#dirtyDot` (7 px, `accent`,
circle). Estados `disabled` por atenuación nativa de icono (no `opacity`).

### D8 — Refresco de tema

`_apply_theme` añade `self.app_bar.refresh_icons(resolved)` (explícito,
antes del loop de toolbars): re-tinta flask/pill/theme/settings/plus/x/doc
con los tokens del tema **resuelto** (light/dark, incluyendo `system`). La
caché del provider lleva el color en la clave → light→dark→light produce
iconos correctos por tema (sin contaminación cruzada). El QSS (applied en
`apply_theme`) cubre fondos/bordes/textos. `action_theme_toggle` queda en
`checked == (resolved == "dark")`.

### D9 — Responsividad

Altura fija 54 px (no crece verticalmente). Pestañas: `elideMode=Right`
(truncado) + `usesScrollButtons=True` (scroll cuando no caben). Pill 250 px
fijos. En 980×600 el ancho mínimo de la app bar (~600 px) cabe con margen.
El canvas (wrapper central) se redimensiona por la VBox; no se rompe el
`QTabWidget`.

## Riesgos y mitigaciones

| Riesgo | Mitigación |
|---|---|
| Drift entre espejo y `QTabWidget` | Resync completo idempotente tras cada cambio estructural (observer) + test de paridad. |
| Autosave dispara al ensuciar en tests | `setClean(False)` → debounce; el test limpia con `setClean(True)` (cancela debounce) al terminar. |
| `QLabel` + `QIcon` (qFatal) | La marca y la pill usan `provider.pixmap(...)`; botones usan `QToolButton`+QAction. |
| Atajo Ctrl+K | No se registra; test aserta que `action_clean_2d_full` conserva su atajo y la pill no tiene shortcut. |
| Orden de toolbars/menús | Menú visible arriba (sin cambios); app bar en central wrapper; `main_toolbar` oculta (no se elimina). |
