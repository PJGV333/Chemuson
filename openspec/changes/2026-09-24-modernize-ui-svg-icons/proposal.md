# Proposal: Sistema de iconos SVG de producción (Fase 2 del PLAN de modernización)

## Why

Todos los iconos de la UI son hoy QPainter programático en
`src/chemuson/gui/icons.py` (~1320 líneas): inconsistentes entre sí (puntero
con relleno, flechas con cabezas a mano, círculos con letras, glifos tipográficos
en Arial), lentos (el tint de undo/redo hace pixel-loop sobre iconos del
sistema) y frágiles al cambio de tema (se regeneran completos en cada
`refresh_icons()`). El spike PyQt6 aprobado (`docs/ui-modernization/pyqt6-spike/`,
commit `59e977d`) y la Fase 1 (tokens + `IconProvider` SVG en
`src/chemuson/gui/theme/icon_provider.py`) ya validaron el sustituto: un set de
SVG 24×24 (trazo `currentColor` 1.75, esquinas redondeadas) teñido por
sustitución de `currentColor` + `QSvgRenderer`, con caché por
(nombre, color, tamaño) y HiDPI. Esta fase convierte ese mecanismo en el
sistema de iconos real de ChemUSON **sin cambiar la estructura de la ventana**
(toolbars, docks, menús y canvas se conservan) y sin romper la API histórica
de `gui/icons.py`.

## What Changes

- **Set de iconos SVG en `src/chemuson/gui/theme/icons/`** (55 archivos
  `i-<name>.svg`): herramientas genéricas (pointer, eraser, pan, rotaciones,
  flip, zoom, chain, lasso, corner, frames, TLC, gel, documentos, undo/redo,
  copy/paste, clean), 11 estilos de enlace, 17 flechas/anotaciones,
  energy-levels, orbital molecular y ancla ondulada. Lenguaje uniforme del
  spike: rejilla 24 px, `viewBox="0 0 24 24"`, trazo `currentColor` 1.75,
  linecap/linejoin redondeados. Set propio del repo (sin descargar
  bibliotecas); atribución y licencia en
  `src/chemuson/gui/theme/icons/LICENSE.txt`.
- **`src/chemuson/gui/theme/icon_svg.py`** (nuevo): generadores de SVG
  dinámico (cadenas, sin Qt) para los iconos parametrizados: átomo CPK
  (texto + color), carga (±), electrones (nº), radical + carga, anillo
  (nº de lados + aromático), plantilla de anillo (etiqueta + lados),
  diagrama de energía (cajas + etiqueta + relleno) y esfera de coordinación
  (color con gradiente radial).
- **`src/chemuson/gui/theme/icon_provider.py`** (extensión aditiva):
  `icon_dynamic(key, color, size, **params)` / `pixmap_dynamic(...)` que
  resuelven el SVG en el registro de `icon_svg` y pasan por el mismo
  render/tinte/caché/HiDPI que el flujo estático. La API de la Fase 1
  (`icon`, `pixmap`, `clear_cache`, `theme_color`, caché por clave) no cambia.
- **`src/chemuson/gui/icons.py` pasa a ser fachada de compatibilidad**:
  mismas firmas públicas (`draw_*`, `get_*_icon`, `set_icon_theme`,
  `icon_*_color`, `ICON_SIZE`, `ATOM_COLORS`), delegando en el
  `IconProvider`. Se elimina el pixel-loop de tint y el uso de
  `QIcon.fromTheme`; los colores por defecto pasan a ser tokens de la Fase 1
  (`icon`, `text3`, `surface3`, `surface`). Los colores químicos/CPK
  (`ATOM_COLORS`, relleno de átomos, esfera de coordinación, pastels de
  niveles de energía) se conservan como colores de dominio.
- **Catálogo `architecture/modules.yml`**: M08 registra la carpeta
  `theme/icons/` y `theme/icon_svg.py`.
- **Tests nuevos** (`tests/test_ui_svg_icons.py`): inventario 1:1 (cada forma/
  enlace/flecha usada por `toolbar.py` y `main_window_ui_builder.py` tiene
  SVG resuelto), validez de viewBox y parseo de cada SVG, provider no nulo en
  16/20/24/28, HiDPI, caché (misma instancia por clave), cambio light→dark→
  light sin iconos con color del tema anterior, glifos dinámicos principales,
  fallo visible sin excepción para asset ausente, compatibilidad de firmas de
  la fachada y smoke de la ventana real en ambos temas.
- **Capturas de referencia** en `docs/ui-modernization/icon-phase-shots/`
  (light/dark de la ventana real post-migración).

## Scope

Fase 2 de `docs/ui-modernization/PLAN.md` ("Sistema de iconos SVG"), con la
restricción explícita de no reorganizar la ventana ni migrar iconos fuera de
la API de `gui/icons.py`.

## Non-goals

- No se crea app bar, document tabs, tool rail nuevo, flyouts, side panel ni
  command palette de producción (Fases 3–6).
- No se migra `gui/orbitals.py` (`draw_orbital_icon`): es un subsistema propio
  con render complejo (lóbulos/orbitales) fuera del alcance de esta fase; se
  deja documentado como residual conocido.
- No se cambian `tool_id`, señales, atajos, handlers, lógica de selección,
  canvas, Clean2D, química, nomenclatura, persistencia ni undo/redo.
- No se cambia el tamaño visual de los iconos en barras existentes (28/24/16
  px según cada toolbar) ni su disposición.
- No se arreglan los fallos preexistentes del baseline (ver `baseline.md`).

## Módulos afectados

- M08 GUI (`src/chemuson/gui/`): `icons.py` (fachada), `theme/` (provider,
  icon_svg, carpeta de assets), `toolbar.py`, `text_toolbar.py`,
  `main_window.py`, `main_window_ui_builder.py` (solo consumers; sin cambios
  de comportamiento).
- `architecture/modules.yml` (registro).
- `tests/`, `docs/ui-modernization/` (evidencia).
