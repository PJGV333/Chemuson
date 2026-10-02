# UI Polish Specification (Fase 7)

## Purpose

Define el pulido final de la UI moderna de ChemUSON (Fase 7 de
`docs/ui-modernization/PLAN.md`): contraste AA del tema claro, estados
deshabilitados inequívocos, tooltips uniformes `Nombre (Shortcut)`, HiDPI
para assets rasterizados, onboarding nativo de primera ejecución, thumbnails
de plantillas escalables, registro de plantillas siempre coherente con la
paleta de comandos y manual actualizado. El canvas sigue siendo hoja blanca;
no se modifica química, geometría, `tool_id`, persistencia ni undo semantics.

## ADDED Requirements

### Requirement: Contraste AA del tema claro

El token `text3` del tema claro SHALL tener contraste mínimo AA (≥ 4.5:1)
sobre las superficies claras donde se usa (fondo de app `#F1F5F9`, surface
`#FFFFFF`). El tema oscuro SHALL conservar sus valores actuales (ya cumplen).

#### Scenario: Token light corregido
- **GIVEN** el tema claro activo
- **WHEN** se resuelven los tokens de diseño
- **THEN** `text3` light es `#5E6E82` (≥ 4.5:1 sobre `#F1F5F9`, `#FFFFFF` y
  `#F8FAFC`)
- **AND** `text2`, `accent`, bordes y todo el bloque dark permanecen
  inalterados.

### Requirement: Estados deshabilitados inequívocos

Los controles deshabilitados del shell SHALL distinguirse visualmente de los habilitados en reposo (opacidad reducida y/o fondo distinto) y el hover SHALL NO hacerlos parecer activos.
Ejemplo: Deshacer sin operaciones deshacibles aparece atenuado y sin resaltado de hover. El estado habilitado SHALL seguir derivado de la `QAction` (`setEnabled`); no se copia ni se fuerza estado en los widgets. Afecta: botones de app bar, botones del rail, celdas de flyout y de la cuadrícula de paletas.

#### Scenario: Deshacer deshabilitado
- **GIVEN** una ventana sin operaciones deshacibles
- **WHEN** el ratón pasa sobre el botón de Deshacer de la app bar
- **THEN** el botón aparece atenuado (opacidad < 1) y sin resaltado de hover
- **AND** el `QAction` de deshacer sigue siendo la fuente de su estado
  `enabled`.

#### Scenario: Celdas de flyout deshabilitadas
- **GIVEN** una celda de flyout o de cuadrícula de paleta deshabilitada
- **WHEN** se pinta la celda
- **THEN** su etiqueta usa `text3`, el fondo usa `surface2` y la opacidad es
  < 1.

### Requirement: Tooltips uniformes `Nombre (Shortcut)`

Los tooltips del rail de herramientas, de la app bar y de las acciones relevantes del panel lateral SHALL seguir la convención `Nombre` o `Nombre (Shortcut)`, usando `QAction.shortcut().toString()` cuando la `QAction` subyacente expone un shortcut no vacío.
En caso contrario se conserva el texto base (incluidos los atajos contextuales de una tecla del rail, que no son `QAction.shortcut()`).

#### Scenario: Tooltip con shortcut de QAction
- **GIVEN** un botón del rail cuya `QAction` histórica tiene shortcut
- **WHEN** se setea el tooltip del botón
- **THEN** el tooltip es `"Nombre (shortcut)"` usando
  `action.shortcut().toString()`.

#### Scenario: Tooltip de Nuevo referenciado a la acción
- **GIVEN** la pestaña de acción "Nuevo" de `DocumentTabBar`
- **WHEN** se setea su tooltip
- **THEN** el tooltip incluye el shortcut de `action_new` (Ctrl+N) en vez de
  una cadena hardcodeada desligada de la acción.

### Requirement: HiDPI sin artefactos en assets rasterizados

Los assets rasterizados de la UI SHALL renderizarse multiplicando por el devicePixelRatio del dispositivo para no aparecer difuminados ni recortados a 125/150/200 % de escala.
Afecta a thumbnails de plantillas y glifos de la barra de formato de texto. El contenido químico del grafo del thumbnail (átomos, enlaces, posiciones relativas) SHALL NO modificarse: solo cambia la escala de render. El `IconProvider` ya es DPR-aware y se conserva.

#### Scenario: Thumbnail a 200 %
- **GIVEN** un entorno offscreen con `QT_SCALE_FACTOR=2`
- **WHEN** se genera el thumbnail de una plantilla con átomos
- **THEN** el `QPixmap` resultante tiene tamaño físico `88×56 × 2` (o mayor)
  con `devicePixelRatio() == 2`
- **AND** el mismo grafo produce el mismo conjunto de átomos y enlaces
  que a `QT_SCALE_FACTOR=1`.

#### Scenario: Smoke HiDPI del shell
- **GIVEN** un entorno offscreen con `QT_SCALE_FACTOR=2`
- **WHEN** se construye la ventana principal
- **THEN** la app bar, el rail de herramientas y el side panel existen y
  están visibles
- **AND** los iconos del rail son `QIcon` no nulos.

### Requirement: Onboarding de primera ejecución

La aplicación SHALL mostrar un overlay de bienvenida nativo Qt **solo la
primera vez** (o tras borrar la preferencia), con exactamente 3 pasos:
(1) Tool Rail — "Elige aquí las herramientas de dibujo y anotación."
(2) Canvas — "Dibuja, selecciona y edita tus estructuras en el lienzo."
(3) SidePanel — "Inspector, validación, propiedades, plantillas y apariencia
están aquí." El overlay SHALL ofrecer `Anterior`, `Siguiente`, `Cerrar` y una
opción "No mostrar de nuevo". La finalización SHALL persistirse mediante
`platform.settings` bajo la clave `ui/onboarding/completed`. El overlay SHALL
ser hijo de la ventana principal (shell) y SHALL NO agregar items a la escena
ni modificar la lógica del canvas.

#### Scenario: Primera ejecución muestra los 3 pasos
- **GIVEN** una `QSettings` sin `ui/onboarding/completed`
- **WHEN** se monta el shell
- **THEN** el overlay se muestra en el paso 1 (Tool Rail)
- **AND** `Siguiente` avanza a Canvas y luego a SidePanel.

#### Scenario: No se repite
- **GIVEN** `ui/onboarding/completed` verdadero
- **WHEN** se monta el shell
- **THEN** el overlay no se muestra.

#### Scenario: "No mostrar de nuevo" persiste
- **GIVEN** el overlay abierto con la opción marcada
- **WHEN** el usuario cierra el onboarding
- **THEN** `ui/onboarding/completed` queda verdadero y el overlay no vuelve a
  mostrarse.

### Requirement: Registro de plantillas coherente en la paleta

Tras una mutación de la biblioteca de plantillas (crear, importar, eliminar o renombrar) la aplicación SHALL refrescar tanto el menú de plantillas como el registro de la paleta de comandos, de modo que la paleta nunca presente plantillas eliminadas ni omita plantillas nuevas cuando se abre.
El mecanismo SHALL reutilizar el refresco existente del menú (fuente de verdad) y la reconstrucción del registro; SHALL NO modificar `TemplateLibrary`, la química de plantillas ni el formato de persistencia.

#### Scenario: Plantilla creada y paleta abierta
- **GIVEN** una ventana con una plantilla recién creada/importada
- **WHEN** se abre la paleta de comandos
- **THEN** la nueva plantilla aparece como entrada ejecutable
- **AND** no aparecen entradas de plantillas eliminadas.

### Requirement: Iconos provisionales de formato de texto sustituidos

La barra de formato de texto SHALL usar iconos SVG del set temático para
alineación (izquierda/centro/derecha) y sub/superíndice, en lugar de glifos
tipográficos. Los botones de negrita/cursiva/subrayado (`B`/`I`/`U`)
conservan su glifo (convención estándar). Los nuevos SVG SHALL integrarse al
`IconProvider` existente (tint por tema, cache, DPR).

#### Scenario: Botones de alineación con SVG
- **GIVEN** la barra de formato de texto visible
- **WHEN** se pintan los botones de alineación y sub/superíndice
- **THEN** sus iconos provienen de `i-align-left`, `i-align-center`,
  `i-align-justify`, `i-subscript` y `i-superscript` (SVG del provider),
  no de glifos tipográficos.

### Requirement: Referencias de atajos actualizadas y manual reflejo de la UI

Los docstrings/comentarios de código activo que documentan la apertura de la paleta de comandos SHALL indicar `Ctrl+P` (la paleta) y `Ctrl+K` únicamente para Clean2D, y el manual SHALL describir la UI moderna con los atajos reales.
El manual (`docs/MANUAL_USUARIO.md`) SHALL describir la UI moderna (AppBar, pestañas de documento, rail de herramientas + flyouts, side panel con tabs, barra de formato de texto, barra de estado, tema claro/oscuro y onboarding) con los atajos reales: `Ctrl+P` (Command Palette), `Ctrl+K` (Limpiar 2D 1 paso), `Ctrl+Shift+K` (Limpiar 2D para publicación), `Ctrl+Alt+K` (Proponer conformero). El manual SHALL NO inventar atajos que no existan.

#### Scenario: Código activo sin stale Ctrl+K de la paleta
- **GIVEN** los módulos `command_palette`, `qss`, `tokens`, `main_window` y
  `app_bar`
- **WHEN** se leen sus docstrings/comentarios de apertura de la paleta
- **THEN** citan `Ctrl+P` para la paleta y `Ctrl+K` solo para Clean2D.

#### Scenario: Manual refleja la UI moderna
- **GIVEN** el manual actualizado
- **WHEN** se consultan las secciones de interfaz
- **THEN** describen AppBar/DocumentTabs/ToolRail/SidePanel/CommandPalette con
  los atajos reales y no describen las barras izquierda/derecha históricas
  como componentes actuales.
