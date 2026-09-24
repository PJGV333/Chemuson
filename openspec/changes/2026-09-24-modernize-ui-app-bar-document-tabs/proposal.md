# Proposal: Barra de aplicación y pestañas de documento (Fase 3 del PLAN de modernización)

## Why

Después de la Fase 1 (tokens + tema) y la Fase 2 (iconos SVG), la aplicación
real funciona y luce moderna en color/iconografía, pero la **estructura de la
ventana sigue siendo la histórica**: QMenuBar + toolbar superior clásica
(iconos new/open/save/undo/redo/rotar/flip/clean/SMI) + `QTabWidget` con su
tab bar nativo sobre el lienzo. El mockup aprobado
(`docs/ui-modernization/mockup-ui.html`) y el spike PyQt6 aprobado
(`docs/ui-modernization/pyqt6-spike/`, commit `59e977d`) definen la primera
transformación estructural visible: una **barra de aplicación** compacta
(54 px) con marca, versión, **pestañas de documento integradas**, botón `+`,
píldora de búsqueda (placeholder de la futura command palette) y controles
derechos (undo/redo, tema, ajustes). Esta fase implementa exactamente ese
shell superior y las pestañas, **reutilizando las QAction reales** (sin
duplizar handlers ni atajos) y el estado real de documentos
(`CanvasTabManager` + `QUndoStack`).

## What Changes

- **`src/chemuson/gui/app_bar.py`** (nuevo): `AppBar` (QFrame, 54 px) con
  marca (flask SVG + "Chemuson" + píldora de versión), `DocumentTabBar`,
  botón `+`, píldora de búsqueda (`SearchPill`, placeholder sin atajo) y
  botones undo/redo/tema/ajustes que **sostienen las QAction existentes**
  (`action_undo`, `action_redo`, `action_preferences`) + una nueva
  `action_theme_toggle` (checkable) conectada al handler existente
  `ChemusonWindow.toggle_theme`.
- **`src/chemuson/gui/document_tabs.py`** (nuevo): `DocumentTabBar` (QTabBar
  con `elideMode=Right`, `drawBase=False`, `expanding=False`, scroll buttons,
  `movable=True`): icono de documento por pestaña, **punto de suciedad**
  (color `accent`, del estado real `undo_stack.isClean()`), botón de cierre
  discreto, botón `+` de esquina, hover/selected con acento inferior
  (lenguaje visual del spike, nativo QtWidgets + QSS de tokens).
- **`src/chemuson/gui/shell/assembly.py`**: el `QTabWidget` y la `AppBar` se
  componen en un wrapper central `QWidget` ( VBox `[app_bar, tabs]`); la tab
  bar nativa del `QTabWidget` se oculta (el widget sigue siendo la fuente de
  verdad de documentos); el `QMenuBar` se **conserva visible y funcional**
  (Archivo/Editar/Ver/Estructura/Reacción/Ayuda + Alt); la toolbar superior
  clásica (`main_toolbar`) se **oculta** (sus acciones siguen accesibles por
  menú, atajo y app bar; el objeto se conserva como compatibilidad).
- **`src/chemuson/gui/tab_manager.py`**: `CanvasTabManager` acepta un
  callback opcional `on_change` (aditivo; por defecto `None`, comportamiento
  idéntico para todos los callers actuales) que se invoca al crear/desechar
  pestañas, para que la `DocumentTabBar` (espejo pasivo) se resincronice.
- **`src/chemuson/gui/main_window.py` + `main_window_ui_builder.py`**:
  sincronización app bar ↔ `QTabWidget` (callbacks `on_change`/`on_tab_updated` del `CanvasTabManager` + señales
  `currentChanged`, `tabMoved`); `action_preferences` se crea en
  `create_local_actions` (mismo objeto, misma conexión, mismo menú);
  `_apply_theme` refresca los iconos de la app bar con el tema resuelto.
- **`src/chemuson/gui/theme/qss.py`**: sección de QSS de tokens para
  `#app_bar`, `#docTabs`, `#tabNewBtn`, `#searchPill`, `#kbdK`,
  `#appBrandName`, `#appVersionPill`, `[appBarBtn]`, `[tabClose]`, `#dirtyDot`.
- **`src/chemuson/gui/theme/icons/`**: 8 SVG nuevos (63 en total): `plus`,
  `search`, `moon`, `sun`, `sliders`, `flask`, `x`, `doc` (geometrías del
  spike aprobado; `currentColor`, 24×24, trazo 1.75).
- **Tests** (`tests/test_ui_app_bar_tabs.py`): app bar/tabs unitarios y de
  ventana real (crear/cerrar/cambiar/reordenar/suciedad/acciones compartidas/
  enabled-disabled/tema light→dark→light/resize 1440×900 y 980×600).
- **Capturas** `docs/ui-modernization/app-bar-phase-shots/` (light/dark,
  1440×900 y 980×600, varias pestañas, pestaña modificada, undo/redo
  enabled/disabled) + script reproducible.
- **Catálogo `architecture/modules.yml`**: M08 registra los nuevos módulos y
  el test.

## Out of Scope (explícito)

- Command palette completa (Fase 6): la píldora de búsqueda es **visual /
  placeholder** y **NO** registra atajo (Ctrl+K sigue siendo "Clean 2D full"
  hoy; reasignarse en su fase).
- Nuevo rail lateral, flyouts, side panel, status bar nueva, canvas nuevo.
- Lógica Clean2D, selección, química, nomenclatura, persistencia molecular,
  reacciones, geometry3d.
- No se eliminan QAction, menús, handlers ni código histórico: la
  `main_toolbar` clásica se oculta (no se borra) y el `QMenuBar` permanece.

## Impact

- Affected specs: `ui-app-bar-document-tabs` (nueva).
- Affected code: M08 (`gui/app_bar.py`, `gui/document_tabs.py`,
  `gui/shell/assembly.py`, `gui/main_window.py`,
  `gui/main_window_ui_builder.py`, `gui/tab_manager.py`), tema
  (`theme/qss.py`, `theme/icons/`), tests, docs.
- Riesgo controlado: el espejo de pestañas es de solo-lectura sobre el
  `QTabWidget` (fuente única de verdad); todo flujo de documento (abrir,
  guardar, cerrar, autosave, recuperación) sigue pasando por
  `CanvasTabManager`/controladores existentes sin cambio de contrato.
