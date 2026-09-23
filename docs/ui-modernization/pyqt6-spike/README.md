# Spike visual PyQt6 — modernización de UI de ChemUSON

Mockup **experimental y descartable** en PyQt6/QtWidgets para evaluar si el lenguaje
visual propuesto en [`../mockup-ui.html`](../mockup-ui.html) y [`../PLAN.md`](../PLAN.md)
es sostenible en Qt **antes** de tocar el código de producción.

> **No es código de producción.** No reutiliza clases internas de ChemUSON ni toca
> `src/chemuson/`, la persistencia, Clean2D, ChemName ni los docks/toolbars reales.
> Todo está aislado en este directorio. Si la decisión es "sí, PyQt6 aguanta el diseño",
> este spike sirve de referencia visual; si no, se elimina sin coste.

## Cómo ejecutarlo

Requiere Python 3 + PyQt6 (en este entorno: PyQt6 6.11, sin `venv`).

```bash
# Ventana interactiva (tema claro, 1440x900)
python docs/ui-modernization/pyqt6-spike/app.py

# Tema oscuro y tamaño a medida
python docs/ui-modernization/pyqt6-spike/app.py --theme dark --size 1440x900

# Modo humo: secuencia de verificación automática + capturas en un directorio
python docs/ui-modernization/pyqt6-spike/app.py --smoke /tmp/spike_shots
```

Opciones:

| Opción | Valores | Por defecto | Descripción |
|---|---|---|---|
| `--theme` | `light`, `dark` | `light` | Tema inicial (conmutable en vivo con el botón de la app bar). |
| `--size` | `ANCHOxALTO` | `1440x900` | Tamaño inicial de la ventana (p. ej. `980x600` para el caso compacto). |
| `--smoke` | `DIR` | — | Ejecuta la secuencia de verificación (temas, flyouts, tabs, paleta, zoom, pestañas) y guarda capturas PNG en `DIR`. Imprime `SMOKE: OK — 0 fallos` si todo pasa. |

En el modo `--smoke` también se puede forzar la plataforma offscreen para CI sin display:

```bash
QT_QPA_PLATFORM=offscreen python docs/ui-modernization/pyqt6-spike/app.py --smoke /tmp/spike_shots
```

## Estructura

| Archivo | Responsabilidad |
|---|---|
| `app.py` | `SpikeWindow` (app bar, tool rail, canvas, panel derecho, status bar), orquestación de la secuencia `--smoke` y `main()`. |
| `theme.py` | Tokens de diseño (color/spacing/radius/typo) para `light`/`dark`, `Theme` (aplica `QPalette` + `QSS` centralizado) y `METRICS` (anchos de paneles, rail, etc.). |
| `widgets.py` | Widgets reutilizables: `ToolRail` (botones hover/pressed/active + flyouts), `Flyout`, `FlyoutCell`, `SideTabRow`/`_SideTabStrip`, `StatusTool`, `SearchPill`, `Kbd`, sombras. |
| `palette.py` | `CommandPalette` (Ctrl+K): overlay hijo con filtro, navegación ↑/↓, Enter/Esc y acciones. |
| `canvas_demo.py` | `QGraphicsView`/`QGraphicsScene` de ejemplo: molécula con rejilla, selección, números y flecha de mecanismo. |
| `icons.py` | `IconProvider`: pipeline SVG → `QIcon`/`QPixmap` (renderiza los `.svg` reales, sin emojis), con `devicePixelRatio` configurable. |
| `icons/*.svg` | Los iconos vectoriales reales del mockup. |

## Diferencias visuales respecto al mockup HTML

Reproducción **buena para decidir**, no pixel-perfect. Diferencias conocidas y honestas:

1. **Sombreado/sombras.** El mockup usa `box-shadow`; aquí se usa `QGraphicsDropShadowEffect`.
   El radio y la difuminación son aproximados y el costó de render es distinto (irrelevante para
   la decisión, pero visible a ojo).
2. **Transparencias offscreen.** En plataforma `offscreen` (CI) Qt no compone la cadena
   `WA_TranslucentBackground` + viewport transparente, así que la tira de tabs del panel
   derecho se pinta con un fondo **opaco** (`surface`) en lugar de translúcido. En un display
   real la transparencia funciona; el color base es el mismo.
3. **Flyouts.** El HTML construye la grilla con CSS `grid`; en Qt se reutiliza un `QGridLayout`
   que **hay que recrear en cada `populate`** (tras vaciarlo, `totalSizeHint()` queda obsoleto en
   Qt 6.11) y la altura se fija a partir de los `sizeHint` de las celdas, porque un layout
   visible+padre no tiene `sizeHint` válido hasta pasar por el event loop. El ancho (244 px),
   el gap (6 px) y el wrap de etiquetas se ajustan a mano para coincidir con el HTML.
4. **Íconos.** Se renderizan los `.svg` reales del mockup vía `QSvgRenderer` → `QPixmap`
   (sin emojis), escalados por `devicePixelRatio`. El tamaño de trazo puede variar un punto
   respecto al `stroke` del CSS.
5. **Doble glifo en el flyout de átomos.** El mockup HTML renderiza el elemento dos veces
   (glifo 13 px bold + etiqueta 10,5 px) porque `FLY_IC` solo define iconos para
   `bond/ring/arrow`. El spike replica ese comportamiento **a propósito** para ser fiel.
6. **Fuente.** El HTML usa la pila del sistema (`-apple-system`, …); Qt usa la fuente
   sans por defecto de la plataforma, que puede no ser idéntica.
7. **Paleta (Ctrl+K).** Centrada con `margin: 11vh auto 0` y `width: min(600, 92vw)`
   como el HTML; implementada como `QWidget` hijo con `eventFilter` (no `QDialog`) para
   que el backdrop cubra la propia app y no aparezca un diálogo del sistema.

## Qué fue fácil en Qt vs. qué requirió widget/QSS propio

**Fácil con `QPalette` + `QSS` (tokens centralizados):**
- Temas claro/oscuro conmutable en vivo (un solo `QSS_TEMPLATE` + `QPalette`).
- App bar, status bar y tool rail con estados hover/pressed/active (selectores `QSS`).
- Tabs de documento, panel derecho con tabs, `SearchPill` y `Kbd` (pill de atajo).
- Paleta de comandos con filtro y navegación.

**Requerir widget/QSS propio:**
- `FlyoutCell` + **wrap de texto manual** (`QLabel` con `wordWrap` no parte palabras más
  anchas que la celda; se implementó `_wrap_text` por espacios/caracteres).
- `_SideTabStrip`: pinta fondo opaco + subrayado de acento en `paintEvent` (la cadena de
  transparencia no se compone en offscreen).
- Sombras vía `QGraphicsDropShadowEffect` (`shadow()`).
- Canvas demo en `QGraphicsView`/`QGraphicsScene` (molécula, rejilla, selección, números).
- `IconProvider` (SVG → pixmap) y el overlay de la paleta (`eventFilter` en vez de `QDialog`).

## Verificación

La secuencia `--smoke` comprueba de forma automatizada: tema claro/oscuro, 2 pestañas de
documento, selector activo, estado del tool, **flyout de enlace (11 ítems) y de átomos
(10 ítems)**, cierre del flyout al elegir, cambio de tabs del panel derecho, filtro de
plantillas, zoom/fit/rejilla/números, abrir y cerrar pestañas, **paleta (27 acciones, filtro
y Esc)**, y que no hay trazas de `Traceback`. Salida esperada: `SMOKE: OK — 0 fallos`.
