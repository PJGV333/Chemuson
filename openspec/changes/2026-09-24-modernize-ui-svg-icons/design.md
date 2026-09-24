# Design: Sistema de iconos SVG de producción

## Contexto

- **Fase 1** ya entregó: tokens (`theme/tokens.py`), QSS/QPalette,
  `theme.apply_theme(...)`, persistencia del tema y `IconProvider`
  (SVG→QIcon/QPixmap, tint por `currentColor`, caché por
  `(name, color, size)`, HiDPI vía `devicePixelRatio`, glifos `glyph:`,
  fallo visible para asset ausente).
- **Callers reales de `gui/icons.py`** (inventario de call sites):
  - `gui/toolbar.py`: `draw_generic_icon` (pointer, lasso, corner, frame,
    rounded_frame, tlc, electrophoresis, chain, rotate_right, undo*, redo*…
    via shapes), `draw_bond_icon` (11 estilos), `draw_arrow_icon` (17 kinds),
    `draw_ring_icon(6, aromatic=True)`, `draw_glyph_icon` (corchetes `[] {}
    () [ { } + - ± δ+ δ- T`, elementos C…H),
    `draw_energy_diagram_icon(slot_count, label_text, label_side, fill_color,
    stroke_visible)`, `draw_energy_levels_icon()`,
    `draw_molecular_orbital_icon()`, `draw_coordination_sphere_icon()`,
    `draw_charge_icon(±)`, `draw_electron_icon(1|2, spread)`,
    `draw_radical_charge_icon(±)`, `draw_wavy_anchor_icon()`.
  - `gui/main_window.py`: `set_icon_theme(resolved_theme)` (en
    `_apply_theme`, antes de `refresh_icons()`).
  - `gui/main_window_ui_builder.py`: `draw_atom_icon("SMI")`,
    `draw_generic_icon` (en `refresh_main_toolbar_icons`).
  - `gui/text_toolbar.py`: `draw_glyph_icon(glyph)` y
    `draw_glyph_icon("■", color=...)`.
- **Tamaños de display**: `ChemusonToolbar`/`SymbolPaletteToolbar` 28 px,
  `main_toolbar` 24 px, `TextToolbar` 16 px (la fachada renderiza a
  `ICON_SIZE = 32` lógicos con dpr; Qt re-escala sin pérdida perceptible).

## Decisiones

### D1. Lenguaje visual: monocromo `currentColor` + opacidades

Todos los iconos de UI usan un único trazo `currentColor` (el tinte lo
proporciona la clave de caché: el token `icon` del tema activo). Los
detalles que antes usaban colores secundarios (`muted`, `fill`, `paper`) se
expresan con `stroke-opacity`/`fill-opacity` sobre `currentColor`, tal como
hizo el spike aprobado. Únicamente se bakan colores literales cuando son
**semánticos de dominio químico**:
- `ATOM_COLORS` (fichas de elemento CPK; el texto blanco/negro se decide por
  luminancia del color, igual que en el código antiguo).
- Esfera de coordinación (color parametrizado + gradiente radial).
- Pastels de los niveles de energía (`#C9EFF8`/`#E6F6AE`/`#F7E6A8`, colores
  de los grupos s/p/d del código actual).
- `fill_color` del diagrama de energía (lo pasa el caller desde el preset).

Esto cumple "lenguaje visual uniforme" del PLAN y mantiene la separación
tokens UI ↔ colores químicos.

### D2. Set propio, sin dependencias

Los 55 SVG son archivos originales del repo escritos a mano en el lenguaje
del spike (que a su vez es compatible con el lenguaje de Lucide/Tabler).
Donde el spike tenía un equivalente aprobado se parte de su geometría
(pointer, undo, redo, chain, lasso, copy, paste, clean, doc-new/open/save,
tlc, gel→electrophoresis, energy→energy-levels, orbital→molecular-orbital,
bond→bond-single, arrow→arrow-forward). `LICENSE.txt` documenta origen y
atribución (set propio; lenguaje inspirado en sets ISC/MIT; sin copiar
archivos de terceros).

### D3. Fachada 1:1 (mismo nombre, misma firma, mismo tipo de retorno)

`icons.py` conserva **exactamente** la API pública actual (ver
`INVENTARIO` más abajo). Cada `draw_*` delega en un `IconProvider`
module-level (lazy, `dpr` desde la app) y devuelve `QIcon`:
- formas/bonds/arrows conocidos → SVG estático `icon(name, color, ICON_SIZE)`;
- parametrizados → `icon_dynamic(key, color, ICON_SIZE, **params)`;
- forma/bond desconocido → `QIcon()` vacío (misma semántica "no dibuja nada"
  del código actual); arrow desconocido → `arrow-line` (el código actual
  cae en una línea simple).
- `set_icon_theme(name)`: `dark`→dark, resto→light (contrato actual) sobre el
  estado `_ICON_THEME`.
- `icon_foreground_color()` → token `icon`; `icon_muted_color()` → token
  `text3`; `icon_fill_color()` → token `surface3`; `icon_paper_color()` →
  token `surface` (del tema activo del icono). Los valores legacy
  (`#111111`/`#F8FAFC`/…) se sustituyen por los tokens de la Fase 1: es el
  cambio visual intencional (tinte slate aprobado en el spike).
- `ICON_SIZE = 32` e `ATOM_COLORS` invariables.

### D4. Dinámicos como SVG generado (sin QPainter)

`icon_svg.py` expone funciones puras `build_<key>(params) -> str` (SVG 24×24)
registradas en un dict. El provider cachea por clave canónica
(`dyn:<key>:<params ordenados>`, color, size). Puntos de fidelidad:
- `ring(sides, aromatic)`: polígono con vértice arriba (ángulo inicial 30°,
  igual que el código QPainter), círculo interno si aromático.
- `atom(text, color)`: círculo relleno CPK + texto negrita; color del texto
  por luminancia (umbral 128, como hoy); 2 letras → fuente menor.
- `energy-boxes(boxes, label, side, fill, stroke)`: 1–N cajas con flecha en
  la caja central; N grande (hasta 56 en presets) produce cajas sub-píxel
  que se fusionan en banda (mismo aspecto que el pixmap actual).
- `sphere(color)`: `<radialGradient>` de 3 paradas (highlight 170 %, base,
  sombra 165 %) + borde.
- `electrons(n, spread)`, `charge(±)`, `radical(±)`: puntos/lineas
  `currentColor` con la misma disposición relativa que hoy.
- `glyph(label)` (ya existe en el provider desde la Fase 1) se usa para
  `draw_glyph_icon`; tamaño de fuente por longitud de etiqueta (1 char →
  mayor), equivalente a la regla actual (15 px/11 px en rejilla 32).

### D5. Semántica de caché y cambio de tema

El color de tinte es **parte de la clave de caché**, así que:
- `set_icon_theme("dark")` + `refresh_icons()` → `draw_*` vuelve a leer el
  color dark → claves nuevas → se renderiza una vez y se cachea; el segundo
  uso devuelve la **misma instancia** (criterio "no redibuja" del PLAN).
- No hay contaminación entre temas: la clave incluye el color, imposible
  servir un icono light en tema dark.
- `clear_cache()` sigue disponible (p. ej. para tests).
- Se elimina `QIcon.fromTheme` + pixel-loop: cero iteración de píxeles.

### D6. Fuera de alcance: `gui/orbitals.py`

`draw_orbital_icon` (M08, ~90 líneas de QPainter por tipo de orbital) no
pertenece a la API de `gui/icons.py`; migrarla exigiría rediseñar lóbulos,
lobos y diagramas específicos. Se documenta como residual conocido; la Fase 4
(flyouts) decidirá su destino. Su import y uso en `toolbar.py` no cambian.

## INVENTARIO 1:1 (API histórica → SVG)

Estáticos (55 archivos `i-<name>.svg` en `theme/icons/`):

| API histórica | SVG |
|---|---|
| `draw_generic_icon("pointer" \| "eraser" \| "pan" \| "lasso" \| "chain" \| "corner" \| "frame" \| "rounded_frame" \| "tlc" \| "electrophoresis" \| "copy" \| "paste" \| "clean")` | `pointer`, `eraser`, `pan`, `lasso`, `chain`, `corner`, `frame`, `rounded-frame`, `tlc`, `electrophoresis`, `copy`, `paste`, `clean` |
| `draw_generic_icon("rotate_left" \| "rotate_right")` | `rotate-left`, `rotate-right` |
| `draw_generic_icon("flip_horizontal" \| "flip_vertical")` | `flip-horizontal`, `flip-vertical` |
| `draw_generic_icon("zoom_in" \| "zoom_out")` | `zoom-in`, `zoom-out` |
| `draw_generic_icon("document_new" \| "document_open" \| "document_save")` | `doc-new`, `doc-open`, `doc-save` |
| `draw_generic_icon("undo" \| "redo")` | `undo`, `redo` (SVG propios; se retira `QIcon.fromTheme`) |
| `draw_bond_icon(<11 estilos>)` | `bond-single`, `bond-bold`, `bond-double`, `bond-triple`, `bond-aromatic`, `bond-interaction`, `bond-coordination`, `bond-wedge`, `bond-hashed`, `bond-wavy`, `bond-flex` |
| `draw_arrow_icon(<17 kinds>)` | `arrow-forward`, `arrow-retro`, `arrow-both`, `arrow-equilibrium`, `arrow-forward-open`, `arrow-retro-open`, `arrow-both-open`, `arrow-equilibrium-open`, `arrow-forward-dashed`, `arrow-retro-dashed`, `arrow-both-dashed`, `arrow-equilibrium-dashed`, `arrow-line`, `arrow-line-dashed`, `arrow-retrosynthetic`, `arrow-curved`, `arrow-curved-fishhook` |
| `draw_energy_levels_icon()` | `energy-levels` (pastels s/p/d bakes) |
| `draw_molecular_orbital_icon()` | `molecular-orbital` |
| `draw_wavy_anchor_icon()` | `wavy-anchor` |
| `draw_generic_icon(<shape desconocido>)` | `arrow-line` (fallback: línea, como hoy) |
| `draw_arrow_icon(<kind desconocido>)` | `arrow-line` (fallback: línea, como hoy) |

Dinámicos (`icon_svg.py`, registro en el provider):

| API histórica | builder(params) |
|---|---|
| `draw_atom_icon(text, color)` | `atom(text, color)` — CPK + texto por luminancia |
| `draw_coordination_sphere_icon(color)` | `sphere(color)` — gradiente radial |
| `draw_glyph_icon(text, color)` | `glyph(label)` (Fase 1) |
| `draw_charge_icon(sign)` | `charge(sign)` |
| `draw_electron_icon(count, spread)` | `electrons(count, spread)` |
| `draw_radical_charge_icon(sign)` | `radical(sign)` |
| `draw_ring_icon(size, aromatic)` | `ring(sides, aromatic)` |
| `draw_ring_template_icon(label, size)` | `ring-template(label, sides)` |
| `draw_energy_diagram_icon(boxes, *, label_text, label_side, fill_color, stroke_visible)` | `energy-boxes(boxes, label, side, fill, stroke)` |

Atajos: `get_pointer_icon/get_eraser_icon/get_single_bond_icon/
get_double_bond_icon/get_benzene_icon` → delegan en `draw_*`.

## Riesgos de regresión y mitigación

1. **`<text>` dentro de SVG**: QSvgRenderer renderiza `<text>` (validado en
   el spike con `i-atom.svg` y en tests de Fase 1 con `glyph:`). Riesgo:
   métrica de fuente entre sistemas. Mitigación: fuente `sans-serif`
   explícita + tests de render no vacío (píxeles no transparentes).
2. **Cambio de color perceptible** (legacy near-black → token slate `icon`):
   intencional (diseño aprobado). Mitigación: capturas light/dark
   comparativas documentadas.
3. **Firma/behavior accidental en la fachada**: el tipo de retorno es
   `QIcon` en todos los casos (como hoy) y las firmas se conservan
   literalmente; test de compatibilidad por `inspect.signature` + llamadas
   reales.
4. **Caché no invalidada al cambiar de tema**: imposible por diseño D5
   (color en la clave); test light→dark→light explícito.
5. **Asset ausente**: `QIcon()`/`QPixmap()` vacíos, sin excepción (contrato
   de Fase 1); el inventario 1:1 + test de paridad sobre el código fuente de
   `toolbar.py` impiden tool_ids huérfanos futuros.
6. **Rendimiento**: el primer render de un (name, color, size) usa
   QSvgRenderer (µs-ms); los menús de paleta construyen sus iconos bajo
   demanda y quedan cacheados; sin pixel-loops.

## Pruebas necesarias (ver `tasks.md`)

- Inventario/paridad: cada literal `draw_generic_icon("x")` /
  `draw_bond_icon("x")` / `draw_arrow_icon("x")` en `toolbar.py` y
  `main_window_ui_builder.py` resuelve a un SVG existente.
- Cada SVG: parse XML válido, `viewBox="0 0 24 24"`, trazo `currentColor`
  (allowlist documentada para colores de dominio bakes).
- Provider: no nulo en 16/20/24/28; HiDPI dpr=2 (ancho físico = 2× size);
  caché = misma instancia; color distinto → pixmap distinto; asset ausente →
  vacío sin excepción.
- Fachada: todas las funciones públicas existen con sus firmas; retorno
  `QIcon` no nulo para todos los nombres del inventario; `set_icon_theme`
  cambia el pixmap; `icon_*_color()` devuelven tokens válidos.
- Dinámicos: anillos 3–8 (aromático/no), átomos C/N/Cl/SMI, cargas ±,
  electrones 1/2, radicales ±, energía 1/3/5/7 con etiqueta y sin, esfera
  con color, glifos `C`/`δ+`/`[]`/`■` (con color).
- Ventana real: creación + `refresh_icons()` en light y dark sin
  excepciones; acciones de toolbar con iconos no nulos.

## Archivos/módulos afectados

- **Nuevos**: `src/chemuson/gui/theme/icons/*.svg` (55) + `LICENSE.txt`,
  `src/chemuson/gui/theme/icon_svg.py`, `tests/test_ui_svg_icons.py`,
  `docs/ui-modernization/icon-phase-shots/*`.
- **Modificados**: `src/chemuson/gui/icons.py` (fachada),
  `src/chemuson/gui/theme/icon_provider.py` (API dinámica aditiva),
  `architecture/modules.yml` (M08).
- **No se modifican**: `toolbar.py`, `text_toolbar.py`, `main_window.py`,
  `main_window_ui_builder.py`, `orbitals.py`, `clean2d/`, `core/`,
  `platform/` (los consumers no cambian: la fachada es 1:1).
