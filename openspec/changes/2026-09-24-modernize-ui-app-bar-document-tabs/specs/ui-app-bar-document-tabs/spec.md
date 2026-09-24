# UI App Bar & Document Tabs Specification

## Purpose

Define la barra de aplicación y las pestañas de documento de Chemuson (Fase
3): un shell superior nativo (QtWidgets + QSS de tokens + iconos SVG del
`IconProvider`) que reutiliza las QAction reales y el estado real de
documentos (`CanvasTabWidget` + `CanvasTabManager` + `QUndoStack`),
aproximando la ventana real al spike aprobado sin introducir la command
palette, el nuevo rail, el side panel ni cambios de química/persistencia.

## ADDED Requirements

### Requirement: Barra de aplicación integrada en la ventana real

La aplicación SHALL exponer una `AppBar` (`src/chemuson/gui/app_bar.py`) de
altura fija 54 px en la parte superior del área central (sobre el
`QTabWidget` de documentos), con: marca (SVG `flask` + nombre "Chemuson" +
píldora de versión), `DocumentTabBar`, botón `+`, píldora de búsqueda
placeholder y botones undo/redo/tema/ajustes a la derecha. El `QMenuBar`
histórico (Archivo, Editar, Ver, Estructura, Reacción, Ayuda) SHALL
permanecer presente y funcional, y la toolbar superior clásica (`main_toolbar`)
SHALL ocultarse sin eliminar sus QAction.

#### Scenario: La ventana real monta la app bar
- **GIVEN** un `ChemusonWindow` instanciado (offscreen)
- **WHEN** se inspecciona la jerarquía de widgets
- **THEN** existe `window.app_bar` (visibilidad `visible`) de altura 54
- **AND** contiene la marca, la `DocumentTabBar`, el botón `+`, la píldora de
  búsqueda y los botones undo/redo/tema/ajustes
- **AND** el widget central es un wrapper que contiene `app_bar` y `tabs`
- **AND** la tab bar nativa de `tabs` está oculta.

#### Scenario: Menú conservado y toolbar clásica oculta
- **GIVEN** el `ChemusonWindow` montado
- **WHEN** se inspeccionan superficies
- **THEN** `menuBar()` es visible con los 6 menús históricos
- **AND** `main_toolbar` existe pero no es visible
- **AND** todas las QAction de `main_toolbar` siguen existiendo en `window`
  y al menos una de (menú, atajo) las mantiene alcanzables.

### Requirement: DocumentTabBar como espejo del estado real de documentos

La `DocumentTabBar` (`src/chemuson/gui/document_tabs.py`) SHALL ser un
espejo de solo-lectura del `QTabWidget` (fuente única de verdad): título por
pestaña desde `CanvasTabManager.tab_titles`, indicador de suciedad desde
`canvas.undo_stack.isClean()`, y sincronización de adición/remoción vía el
observer `on_change` de `CanvasTabManager`. Abrir, cerrar, reordenar y
seleccionar documentos SHALL seguir utilizando el flujo existente
(controladores + `CanvasTabManager` + QAction), sin un segundo estado de
suciedad.

#### Scenario: Pestaña inicial
- **GIVEN** el `ChemusonWindow` recién creado
- **WHEN** se inspecciona la `DocumentTabBar`
- **THEN** contiene exactamente 1 pestaña
- **AND** su texto es "Sin título"
- **AND** su punto de suciedad está oculto
- **AND** está seleccionada.

#### Scenario: Documento nuevo por el botón +
- **GIVEN** el `ChemusonWindow` con 1 pestaña
- **WHEN** se pulsa el botón `+` de la app bar
- **THEN** `window.action_new` se dispara (mismo QAction)
- **AND** el `QTabWidget` y la `DocumentTabBar` quedan con 2 pestañas
- **AND** la pestaña nueva es la activa (canvas nuevo activo).

#### Scenario: Indicador de suciedad desde el estado real
- **GIVEN** una pestaña limpia
- **WHEN** se ejecuta `canvas.undo_stack.setClean(False)`
- **THEN** el `QTabWidget` añade el sufijo `" *"` (contrato existente de
  `update_tab_title`)
- **AND** la `DocumentTabBar` muestra su punto de suciedad
- **AND** cuando se ejecuta `setClean(True)` el punto se oculta de nuevo.

#### Scenario: Cerrar documento reutiliza el flujo existente
- **GIVEN** dos pestañas limpias
- **WHEN** se pulsa el botón de cierre de la pestaña activa
- **THEN** se invoca el flujo `_on_tab_close_requested`/`_close_canvas_tab`
  existente (sin diálogo de confirmación al estar limpias)
- **AND** quedan 1 pestaña en el `QTabWidget` y en la `DocumentTabBar`.

#### Scenario: Cambiar y reordenar pestañas
- **GIVEN** dos pestañas
- **WHEN** se pulsa la pestaña 2 del espejo
- **THEN** `window.canvas` pasa al canvas de la pestaña 2 (handler existente)
- **WHEN** se emite `tabMoved(0, 1)` en el espejo
- **THEN** el `QTabWidget` reordena sus pestañas (contrato `moveTab`)
  y la selección se mantiene coherente.

### Requirement: Reutilización de QAction y estado habilitado

Los botones de la app bar SHALL sostener las QAction existentes
(`action_undo`, `action_redo`, `action_preferences`) mediante
`setDefaultAction` (mismo objeto, sin duplicar handlers ni atajos) y
reflejar su estado habilitado; el botón de tema SHALL usar el handler
existente `toggle_theme`. La píldora de búsqueda SHALL ser un placeholder
sin atajo registrado (Ctrl+K sigue perteneciendo a
`action_clean_2d_full`).

#### Scenario: Mismos QAction compartidos
- **GIVEN** el `ChemusonWindow` montado
- **WHEN** se inspeccionan los botones de la app bar
- **THEN** el botón undo tiene `defaultAction() is window.action_undo`
- **AND** el botón redo tiene `defaultAction() is window.action_redo`
- **AND** el botón ajustes tiene `defaultAction() is window.action_preferences`
- **AND** `action_undo`/`action_redo` conservan sus atajos originales
  (atajos estándar Undo/Redo de la plataforma) sin duplicados.

#### Scenario: Estado enabled real de undo/redo
- **GIVEN** un documento vacío (sin undo disponible)
- **WHEN** se inspecciona el botón undo de la app bar
- **THEN** está deshabilitado (igual que `action_undo`)
- **WHEN** `action_undo.setEnabled(True)`
- **THEN** el botón pasa a habilitado sin lógica adicional.

#### Scenario: Placeholder de búsqueda sin atajo
- **GIVEN** la app bar montada
- **WHEN** se inspecciona la píldora de búsqueda
- **THEN** no registra ningún `QShortcut`/atajo
- **AND** `action_clean_2d_full` conserva su atajo Ctrl+K.

### Requirement: Tema light/dark en app bar y pestañas

La app bar, sus iconos, la `DocumentTabBar` (hover, selected, dirty, close)
y la píldora SHALL actualizarse al cambiar de tema light → dark → light,
sin que ningún icono conserve el tinte del tema anterior (caché del
`IconProvider` con el color en la clave) y sin crecer verticalmente.

#### Scenario: Ciclo de tema sin contaminación
- **GIVEN** el `ChemusonWindow` en tema claro
- **WHEN** se cambia a oscuro (`toggle_theme`)
- **THEN** los iconos de la app bar se re-tinean con tokens oscuros
  (p. ej. la marca usa `accent` oscuro y los botones `icon` oscuro)
- **AND** el botón de tema muestra `sun`
- **WHEN** se vuelve a claro
- **THEN** todos los iconos vuelven a los tintes claros exactos
- **AND** el botón de tema muestra `moon`
- **AND** la altura de la app bar permanece 54 px.

### Requirement: Responsividad mínima

Con 1440×900 y 980×600 la ventana SHALL renderizar sin excepciones: la app
bar mantiene su altura, las pestañas se truncan (elide) o desplazan
(scroll buttons) según el espacio, los controles esenciales (undo/redo,
tema, ajustes, `+`) permanecen presentes y el área del canvas sigue
funcional.

#### Scenario: Redimensionado a tamaños objetivo
- **GIVEN** el `ChemusonWindow` con varias pestañas
- **WHEN** se redimensiona a 1440×900 y luego a 980×600 (offscreen)
- **THEN** no se lanzan excepciones
- **AND** la app bar mantiene 54 px de altura en ambos tamaños
- **AND** la `DocumentTabBar` sigue presente con sus pestañas
- **AND** el `QTabWidget` sigue siendo el área de documentos central.

### Requirement: Iconos SVG de la app bar válidos

Los 8 SVG nuevos (`plus`, `search`, `moon`, `sun`, `sliders`, `flask`, `x`,
`doc`) SHALL cumplir el contrato del set de la Fase 2 (XML válido,
`viewBox="0 0 24 24"`, `currentColor`, sin raster) y estar cubiertos por el
inventario/paridad de `tests/test_ui_svg_icons.py`.

#### Scenario: Inventario ampliado
- **GIVEN** la carpeta `src/chemuson/gui/theme/icons/`
- **WHEN** se valida el inventario completo (63 archivos)
- **THEN** los 8 nombres nuevos existen y son válidos
- **AND** `IconProvider` resuelve cada uno como `QIcon` no vacío en
  16/20/24/28 px.
