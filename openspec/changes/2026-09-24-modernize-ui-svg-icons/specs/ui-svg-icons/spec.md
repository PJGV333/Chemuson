# UI SVG Icons Specification

## Purpose

Define el sistema de iconos de la UI de Chemuson sobre un set de SVG
uniformes (24×24, trazo `currentColor`) servidos por un `IconProvider`
theme-aware, HiDPI y con caché, manteniendo la compatibilidad total con la
API histórica de `gui/icons.py`. Esta especificación no cambia la estructura
de la ventana, los `tool_id`, las señales, los atajos, la lógica de
selección, el canvas, Clean2D, la química, la persistencia ni el undo/redo.

## ADDED Requirements

### Requirement: Set de iconos SVG uniforme y versionado en el repo

El sistema de iconos SHALL almacenar sus iconos estáticos como archivos SVG
en `src/chemuson/gui/theme/icons/` con nombre `i-<name>.svg`, cada uno con
`viewBox="0 0 24 24"`, trazo principal `currentColor` de grosor 1.75 y
esquinas/tamaños redondeados, sin imágenes rasterizadas y sin dependencias
de runtime externas. Los colores literales no `currentColor` solo se
permiten cuando son semánticos de dominio químico (CPK, niveles de energía)
y deberán estar documentados.

#### Scenario: Inventario completo y válido
- **GIVEN** la carpeta `src/chemuson/gui/theme/icons/`
- **WHEN** se itera sobre cada archivo `i-*.svg`
- **THEN** cada archivo es XML/SVG válido
- **AND** declara `viewBox="0 0 24 24"`
- **AND** usa `currentColor` para su trazo principal
- **AND** no hay archivos `.png` ni otros formatos raster en la carpeta.

#### Scenario: Paridad 1:1 con la API histórica
- **GIVEN** el inventario de la Fase 2 (formas genéricas, 11 enlaces,
  17 flechas, energía, orbital, ancla)
- **WHEN** cada nombre del inventario se resuelve contra la carpeta
- **THEN** existe un `i-<name>.svg` para cada nombre
- **AND** no existe ningún `tool_id`/forma usado por `toolbar.py` o
  `main_window_ui_builder.py` que quede sin SVG.

### Requirement: IconProvider dinámico aditivo

El `IconProvider` SHALL ofrecer, además de la API estática de la Fase 1
(`icon`, `pixmap`, `clear_cache`, `theme_color`), una API dinámica
`icon_dynamic(key, color, size, **params)` y `pixmap_dynamic(...)` que
resuelve el SVG a partir de un registro de generadores puros
(`icon_svg.BUILDERS`), lo tiñe con `color` y lo cachea con una clave que
incluye el `key` y los `params` canónicos, el `color` y el `size`. El
comportamiento de la API estática SHALL mantenerse idéntico.

#### Scenario: Icono dinámico parametrizado
- **GIVEN** un `IconProvider` y un builder registrado (p. ej. `ring`)
- **WHEN** se pide `icon_dynamic("ring", color, size, sides=6, aromatic=True)`
- **THEN** se devuelve un `QIcon` no nulo
- **AND** repetir la misma llamada devuelve la misma instancia cacheada
- **AND** cambiar `sides` o `aromatic` produce un icono distinto.

#### Scenario: API estática inalterada
- **GIVEN** un `IconProvider`
- **WHEN** se usan `icon`, `pixmap`, `clear_cache` y `theme_color`
- **THEN** devuelven los mismos tipos y semántica que la Fase 1
  (QIcon/QPixmap, caché por `(name, color, size)`, fallo visible en ausente).

### Requirement: Fachada de compatibilidad en `gui/icons.py`

El módulo `src/chemuson/gui/icons.py` SHALL mantener la misma API pública
histórica (mismos nombres, firmas y tipo de retorno `QIcon`): `set_icon_theme`,
`icon_foreground_color`, `icon_muted_color`, `icon_fill_color`,
`icon_paper_color`, `ICON_SIZE`, `ATOM_COLORS` y las funciones `draw_*`
(generic, bond, arrow, ring, ring_template, atom, glyph, charge, electron,
radical_charge, energy_diagram, energy_levels, molecular_orbital,
coordination_sphere, wavy_anchor) y `get_*_icon`. Cada función SHALL delegar
en el `IconProvider` y SHALL dejar de pintar con QPainter manual, iterar
píxeles o usar `QIcon.fromTheme`. Los colores por defecto de los iconos
monocromos SHALL derivar de los tokens de la Fase 1; los colores CPK de
`ATOM_COLORS` y los rellenos químicos SHALL conservarse como colores de
dominio.

#### Scenario: Mismo contrato de firmas
- **GIVEN** la fachada reescrita
- **WHEN** se inspeccionan las firmas públicas
- **THEN** cada nombre histórico existe
- **AND** acepta los mismos argumentos posicionales/por defecto que la
  versión anterior
- **AND** devuelve `QIcon`.

#### Scenario: Iconos no nulos para el inventario completo
- **GIVEN** una `QApplication`
- **WHEN** se llama a cada `draw_*` con cada valor del inventario
  (24 formas genéricas, 11 enlaces, 17 flechas, anillos 3–8, átomos
  C/N/O/S/P/F/Cl/Br/SMI, cargas ±, electrones 1/2, radicales ±)
- **THEN** cada llamada devuelve un `QIcon` no nulo
- **AND** no se lanza ninguna excepción.

#### Scenario: Formas desconocidas fallan de forma segura
- **GIVEN** la fachada
- **WHEN** se llama `draw_generic_icon("<shape inexistente>")` o
  `draw_bond_icon("<tipo inexistente>")`
- **THEN** se devuelve un `QIcon` no nulo en blanco sin excepción
- **AND** cuando se llama `draw_arrow_icon("<kind inexistente>")` se devuelve
  un `QIcon` no nulo con la línea simple (mismo fallback de la versión
  anterior).

#### Scenario: Sin iteración de píxeles ni iconos de sistema
- **GIVEN** la fachada y el provider
- **WHEN** se generan los iconos de undo/redo y cualquier icono teñido
- **THEN** no se usa `QIcon.fromTheme`
- **AND** no se itera sobre píxeles para recolorear.

### Requirement: Tinte theme-aware sin pixel-loop

El sistema SHALL teñir los iconos sustituyendo `currentColor` en el SVG y
renderizando con `QSvgRenderer` (sin iterar píxeles). El color de tinte de
un icono monocromo SHALL ser parte de la clave de caché, de modo que el
cambio de tema produzca iconos correctos por tema y un uso repetido dentro
del mismo tema no re-renderice.

#### Scenario: Cambio light → dark → light
- **GIVEN** un icono monocromo (p. ej. pointer) bajo tema light
- **WHEN** se conmuta a `set_icon_theme("dark")` y se regenera el icono
- **THEN** el pixmap resultante usa el color dark
- **AND** no contiene el color light anterior (sin caché contaminada)
- **AND** volver a light reproduce el pixmap light original.

#### Scenario: Sin iteración de píxeles ni iconos de sistema
- **GIVEN** la fachada y el provider
- **WHEN** se generan los iconos de undo/redo y cualquier icono teñido
- **THEN** no se usa `QIcon.fromTheme`
- **AND** no se itera sobre píxeles para recolorear.

### Requirement: HiDPI y tamaños estándar

El provider SHALL renderizar con el `devicePixelRatio` de la aplicación y
SHALL producir pixmaps no nulos para los tamaños estándar 16, 20, 24 y 28.

#### Scenario: HiDPI
- **GIVEN** un `IconProvider` construido con `dpr=2`
- **WHEN** se pide un pixmap de tamaño 24
- **THEN** el pixmap físico tiene ancho 48 y `devicePixelRatio() == 2`.

#### Scenario: Tamaños estándar no nulos
- **GIVEN** un icono estático del set
- **WHEN** se solicita en 16, 20, 24 y 28
- **THEN** cada uno devuelve un `QIcon`/`QPixmap` no nulo.

## OUT OF SCOPE

- `gui/orbitals.py` (`draw_orbital_icon`): residual QPainter conocido, fuera
  de esta fase (ver design.md D6).
- Estructura de la ventana, `tool_id`, señales, atajos, canvas, Clean2D,
  química, persistencia, undo/redo: no cambian.
