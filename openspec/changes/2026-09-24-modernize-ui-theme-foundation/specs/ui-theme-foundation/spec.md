# UI Theme Foundation Specification

## Purpose

Define la fundación visual de la UI de Chemuson: un sistema centralizado de
design tokens (claro/oscuro), la generación de QSS/QPalette a partir de esos
tokens, el punto central de aplicación del tema, la persistencia mínima de la
elección de tema y la infraestructura de iconos SVG preparada para la fase de
migración de iconos. Esta especificación no cambia la estructura de la
ventana ni la lógica del canvas.

## ADDED Requirements

### Requirement: Design tokens centralizados por tema

El sistema de temas SHALL definir una tabla de tokens para los temas `light`
y `dark` que cubra, como mínimo: fondo de aplicación (`bg`), superficie
(`surface`), superficie secundaria (`surface2`), superficie terciaria
(`surface3`), bordes (`border`, `borderStrong`), jerarquía de texto
(`text1`, `text2`, `text3`), acento (`accent`, `accentStrong`, `accentHover`,
`accentSoft`, `accentBorder`, `onAccent`), severidades (`danger`, `dangerSoft`,
`warn`, `warnSoft`, `ok`, `okSoft`) y los tokens de dominio `sheet`,
`sheetGrid` y `canvasBg`. Los valores del tema aprobado del spike
(`docs/ui-modernization/pyqt6-spike/theme.py`, commit `59e977d`) SHALL ser la
referencia de estos valores.

#### Scenario: Resolución de tokens light/dark
- **GIVEN** el sistema de tokens está importado
- **WHEN** se resuelven los tokens para `light` y para `dark`
- **THEN** cada tema expone todos los tokens requeridos
- **AND** `bg` de `light` es `#F1F5F9` y de `dark` es `#0B1120`
- **AND** `accent` de `light` es `#0E7490` y de `dark` es `#22D3EE`.

#### Scenario: Colores de tokens válidos
- **GIVEN** la tabla de tokens de cualquiera de los dos temas
- **WHEN** cada valor de color se convierte a `QColor`
- **THEN** el `QColor` resultante es válido (`isValid()`) para todos los tokens de color.

### Requirement: Métricas UI centralizadas

El sistema de temas SHALL centralizar las métricas de UI: espaciado (grilla
de 8 px con escalones 4/8/12/16), radios (superficies 10, botones 9, chips 8)
y tipografía (base 13 px, secundaria 12 px, pequeña 11 px), de forma que el
QSS y los widgets futuros las consuman desde el mismo módulo de tokens.

#### Scenario: Métricas presentes
- **GIVEN** el módulo de tokens
- **WHEN** se consultan las métricas
- **THEN** los espaciados, radios y tamaños de fuente definidos por este
  requirement están disponibles y son numéricos.

### Requirement: Hojas de estilo generadas desde tokens

La UI SHALL generar sus hojas de estilo (`get_main_stylesheet(theme_name)` y
`get_tool_palette_stylesheet(theme_name)`) exclusivamente a partir de la
tabla de tokens del tema, sin colores de UI hardcoded fuera de tokens. Los
colores químicos de dominio (CPK, hoja del canvas) pueden seguir viviendo en
los subsistemas que ya los gobiernan (`gui/icons.py`, canvas), que no se
modifican en esta fase.

#### Scenario: QSS consume tokens
- **GIVEN** los generadores de QSS del tema `dark`
- **WHEN** se generan las hojas de estilo
- **THEN** los valores de `bg`, `surface` y `accent` del tema dark aparecen
  en el QSS generado
- **AND** ningún valor hexadecimal de UI ajeno a la tabla de tokens es
  interpolado como color de interfaz.

#### Scenario: Estados de widget coherentes
- **GIVEN** el QSS principal de cualquiera de los dos temas
- **WHEN** se inspeccionan las reglas de `QToolButton`, `QPushButton`,
  `QLineEdit` y `QTabBar`
- **THEN** los estados `hover`, `pressed`/`focus` y `checked`/`selected`
  usan tokens de superficie y acento (`surface2`, `surface3`, `accentSoft`,
  `accentBorder`, `accentHover`)
- **AND** los estados `disabled` usan tokens de texto/borde atenuados
  (`text3`, `border`, `borderStrong`).

### Requirement: Punto central de aplicación del tema

La aplicación del tema SHALL hacerse desde una única API central
`apply_theme(target, theme_name)` que aplique, en este orden, la fuente base
del tema, el `QPalette` construido desde tokens y la hoja de estilo principal
a una `QApplication` o a la ventana principal. La API SHALL aceptar los
nombres `light`, `dark` y `system`; `system` SHALL resolverse al esquema del
sistema operativo mediante `QStyleHints` cuando esté disponible, y cualquier
nombre desconocido SHALL normalizarse a `light` sin lanzar excepciones.

#### Scenario: Aplicación de ambos temas sin excepciones
- **GIVEN** una `QApplication` (offscreen) y la API central de temas
- **WHEN** se aplica `light` y después `dark`
- **THEN** ninguna excepción se lanza
- **AND** el stylesheet resultante contiene los valores de tokens del tema
  aplicado
- **AND** la paleta del destino coincide con los tokens (p. ej. `Window` =
  `bg` del tema).

#### Scenario: Nombre desconocido no rompe
- **GIVEN** la API central de temas
- **WHEN** se aplica un nombre de tema desconocido (`"neon"`)
- **THEN** el sistema resuelve `light` sin excepción.

#### Scenario: Preparación para "seguir sistema"
- **GIVEN** la API central de temas
- **WHEN** se aplica `"system"` o se invoca `set_theme_from_system(target)`
- **THEN** el tema aplicado es `light` o `dark` (nunca `system` como valor
  final) y coincide con el esquema del sistema cuando `QStyleHints`
  está disponible.

### Requirement: Compatibilidad de la API de estilos existente

El módulo `chemuson.gui.styles` SHALL seguir exportando `get_main_stylesheet`,
`get_tool_palette_stylesheet`, `DEFAULT_THEME`, `MAIN_STYLESHEET`,
`TOOL_PALETTE_STYLESHEET`, `LIGHT_COLORS` y `DARK_COLORS`, de modo que los
callers existentes (`main_window.py`, `toolbar.py`) no requieran cambios.
Los generadores de la fachada SHALL ser equivalentes a los del nuevo sistema
de temas.

#### Scenario: Fachada compatible
- **GIVEN** el módulo `chemuson.gui.styles` importado
- **WHEN** se consultan sus nombres públicos
- **THEN** `get_main_stylesheet("light")` y `get_main_stylesheet("dark")`
  devuelven hojas de estilo no vacías idénticas a las del nuevo sistema
- **AND** `MAIN_STYLESHEET` y `TOOL_PALETTE_STYLESHEET` son strings no vacíos
- **AND** `LIGHT_COLORS` y `DARK_COLORS` contienen las claves legadas
  existentes.

### Requirement: Persistencia mínima de la elección de tema

La preferencia de tema SHALL persistirse mediante el módulo de plataforma
(M21) bajo la clave `ui/theme` con valores `light`, `dark` o `system`;
cualquier valor ausente o inválido SHALL resolverse a `light`. La carga al
arranque SHALL reemplazar el valor hardcodeado actual, y la elección en
Preferencias SHALL persistirse al aplicarse. M21 SHALL continuar sin importar
widgets de GUI.

#### Scenario: Round-trip de preferencias UI
- **GIVEN** un store de preferencias con `ui/theme` = `"dark"`
- **WHEN** se cargan las preferencias UI
- **THEN** el tema resultante es `"dark"`
- **AND** al guardar preferencias con `theme = "system"` la clave `ui/theme`
  queda `"system"`.

#### Scenario: Valor inválido normalizado
- **GIVEN** un store con `ui/theme` = `"neon"`
- **WHEN** se cargan las preferencias UI
- **THEN** el tema resultante es `"light"`.

#### Scenario: Arranque carga el tema persistido
- **GIVEN** una preferencia persistida `ui/theme` = `"dark"`
- **WHEN** se ensambla la ventana principal
- **THEN** `current_theme` es `"dark"` y el tema oscuro se aplica.

### Requirement: Infraestructura de iconos SVG preparada

El paquete de temas SHALL incluir un `IconProvider` que convierta SVG en
`QIcon`/`QPixmap` con tinte por sustitución de `currentColor` (sin iterar
píxeles), caché por `(nombre, color, tamaño)` —de modo que el mismo `QIcon`
se reutilice al re-consultar con la misma clave, lo que hace el provider
theme-aware— y soporte HiDPI vía `devicePixelRatio`. Los iconos SVG
estáticos SHALL buscarse en `src/chemuson/gui/theme/icons/` y los nombres
ausentes SHALL devolver iconos/pixmaps vacíos sin excepciones. La migración
de los iconos de `gui/icons.py` a este provider queda explícitamente fuera
de esta fase.

#### Scenario: Caché por clave
- **GIVEN** un `IconProvider`
- **WHEN** se solicita el mismo icono (nombre, color y tamaño) dos veces
- **THEN** se devuelve la misma instancia de `QIcon` (hit de caché).

#### Scenario: Render HiDPI
- **GIVEN** un `IconProvider` con `devicePixelRatio` 2.0
- **WHEN** se genera un pixmap de 20 px
- **THEN** el pixmap físico mide 40 px y reporta `devicePixelRatio` 2.0.

#### Scenario: Fallo visible sin excepción
- **GIVEN** un `IconProvider`
- **WHEN** se solicita un nombre de icono inexistente
- **THEN** se devuelve un `QIcon` vacío y un `QPixmap` vacío sin lanzar
  excepciones.

### Requirement: Sin regresión estructural

Esta fase SHALL no modificar la lógica de `clean2d/`, `chemname/`,
`chemio/persistence.py`, el canvas (`gui/canvas/`, M09), los comandos
undo/redo, la selección, la geometría, las reacciones ni el formato `.cmsn`,
y SHALL no reorganizar la ventana (menús, toolbars, docks y tabs actuales se
mantienen; solo reciben el nuevo QSS).

#### Scenario: Suite sin regresiones
- **GIVEN** la baseline de la suite completa registrada antes del cambio
- **WHEN** se ejecuta la suite completa tras el cambio
- **THEN** no hay nuevas fallas contra la baseline
- **AND** `git diff` no toca archivos de `clean2d/`, `chemname/`,
  `chemio/persistence.py`, `gui/canvas/`, `gui/icons.py`, `gui/toolbar.py`,
  `gui/docks.py` ni `gui/style.py`.
