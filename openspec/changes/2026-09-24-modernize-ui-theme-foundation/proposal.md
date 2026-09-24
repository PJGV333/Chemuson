# Proposal: UI Theme Foundation and Design Tokens (Fases 0–1 del PLAN de modernización)

## Why

El look de la aplicación no está gobernado por design tokens: los colores viven
dispersos en `src/chemuson/gui/styles.py` (paletas plano QSS), `src/chemuson/gui/icons.py`
(colores QPainter propios) y no hay métricas de espaciado, radios ni tipografía
centralizadas. El spike PyQt6 aprobado (`docs/ui-modernization/pyqt6-spike/`,
commit `59e977d`, `SMOKE: OK — 0 fallos`) demostró que una tabla de tokens
genera un QSS/QPalette visualmente equivalente al mockup. Esta fase convierte
ese resultado en la **fuente de verdad visual de producción** sin cambiar la
estructura de la ventana (app bar, rail, side panel y command palette son
fases posteriores).

## What Changes

- **Nuevo subpaquete `src/chemuson/gui/theme/`** (dentro de M08):
  - `tokens.py`: tabla de tokens light/dark (colores + métricas de espaciado,
    radios, tipografía) extraída del spike aprobado.
  - `qss.py`: generadores `get_main_stylesheet(theme)` y
    `get_tool_palette_stylesheet(theme)` construidos 100 % desde tokens.
  - `palette.py`: `build_qpalette(theme)` (QPalette desde tokens) y
    `theme_font(...)` (fuente base 13 px con fallback de familias).
  - `icon_provider.py`: `IconProvider` SVG→QIcon/QPixmap theme-aware, HiDPI y
    con caché (infraestructura mínima para la Fase 2; sin migrar aún los
    iconos de `icons.py`).
  - `__init__.py`: API central `apply_theme(target, theme_name)`,
    `resolve_theme_name(name)` (prepara "seguir sistema" vía `QStyleHints`),
    `set_theme_from_system(target)`.
- **`styles.py` pasa a ser fachada de compatibilidad**: re-exporta
  `get_main_stylesheet`, `get_tool_palette_stylesheet`, `DEFAULT_THEME`,
  `MAIN_STYLESHEET`, `TOOL_PALETTE_STYLESHEET`, `LIGHT_COLORS`, `DARK_COLORS`
  (los dos últimos como alias deprecation hacia tokens; no hay callers en el
  repo fuera de `styles.py`).
- **Persistencia del tema (M21)**: `UiPreferences` + `load_ui_preferences` /
  `save_ui_preferences` (clave `ui/theme`) en `platform/settings.py`, export
  en `platform/__init__.py`. M21 sigue sin importar QtWidgets/GUI.
- **Aplicación central del tema**: `ChemusonWindow._apply_theme()` delega en
  `theme.apply_theme(...)` (QSS + QPalette + fuente) y conserva el QSS
  específico de toolbars y el refresco de iconos; `shell/assembly.py` carga el
  tema persistido en lugar de hardcodear `"light"`; `_apply_preferences`
  persiste el tema elegido.
- **Catálogo `architecture/modules.yml`**: M08 registra `src/chemuson/gui/theme/`
  (paths + `internal_api: theme` + test nuevo); M21 registra la nueva API
  pública de preferencias UI.
- **Tests nuevos**: resolución de tokens light/dark, validez de colores,
  aplicación de ambos temas sin excepciones, smoke de ventana real,
  compatibilidad de la fachada `styles.py`, caché/HiDPI del IconProvider y
  round-trip de preferencias UI en M21.

## Scope

Fase 0 (preparación: OpenSpec, baseline, registro arquitectónico) + Fase 1
(tokens y QSS) de `docs/ui-modernization/PLAN.md`, con la restricción
explícita de no reorganizar la ventana.

## Non-goals

- No se implementa app bar, pestañas de documento unificadas, rail, flyouts,
  side panel en tabs, barra de estado nueva ni command palette (Fases 3–6).
- No se migran los iconos de `icons.py` a SVG (Fase 2): `icons.py` y
  `set_icon_theme` se mantienen intactos; `IconProvider` se entrega como
  infraestructura preparada y testeada.
- No se añade "Seguir sistema" a la UI (combo de Preferencias): la API
  (`resolve_theme_name("system")`, persistencia `ui/theme=system`) queda
  preparada y documentada; la opción visible llega con la Fase de
  preferencias/ajustes.
- Sin cambios en `clean2d/`, `chemname/`, `chemio/persistence.py`, canvas
  (`gui/canvas/`, M09), comandos undo/redo, selección, geometría, reacciones,
  `.cmsn` ni menús/toolbars/docks existentes (solo reciben el nuevo QSS).
- Sin nuevas dependencias externas (PyQt6/QtSvg ya están en el entorno).

## Compatibility

- `chemuson.gui.styles` conserva sus nombres públicos actuales (los dos
  generadores, `DEFAULT_THEME`, `MAIN_STYLESHEET`, `TOOL_PALETTE_STYLESHEET`
  y las dos paletas legadas como alias de tokens).
- `tool_id`, señales `tool_changed(str)`, jerarquía de mixins de
  `main_window.py` y orden de construcción del shell: intactos.
- M21: extensión puramente aditiva (`UiPreferences` y load/save); `QSettings`
  legacy sin la clave `ui/theme` resuelve a `"light"` (comportamiento
  actual).

## Success criteria

- `openspec validate 2026-09-24-modernize-ui-theme-foundation --strict` pasa.
- Los dos temas se aplican desde un único punto central sin excepciones;
  light/dark visibles en la ventana real (capturas offscreen en
  `docs/ui-modernization/foundation-shots/`).
- Ningún color de UI hardcoded fuera del sistema de tokens salvo excepciones
  químicas justificadas (CPK/canvas, en `icons.py` y canvas: intactos esta
  fase).
- Suite completa sin regresiones contra la baseline; `ruff` scoped limpio;
  `git diff --check` limpio; commit único en `ui/modernization`.

## Discrepancias encontradas en la auditoría (registro)

1. **PLAN.md §1.2 está desactualizado**: el mapa de archivos de `gui/` omite
   módulos actuales (`template_*.py`, `energy_diagrams.py`, `orbitals.py`,
   `periodic_table.py`, `plate_items.py`, `semantic_diagram_workflow.py`,
   `items.py`, `wedge_geometry.py`, `geom.py`, ...). No bloquea; se usa el
   código real como fuente de verdad.
2. **Tabla de tokens de PLAN.md §2.2 vs spike aprobado**: el spike (referencia
   visual aprobada) usa vocabulario camelCase (`bg`, `surface`, `borderStrong`,
   `text1`...) y valores que difieren ligeramente de la tabla de PLAN (p. ej.
   dark `surface2` `#16213A` vs `#1E293B`; dark `text2` `#C3CEDF` vs
   `#CBD5E1`; PLAN omite `accentHover`/`accentBorder`/`onAccent`/`canvasBg`).
   Decisión: **el spike es la referencia** (validado visualmente); la tabla de
   tokens de producción copia el spike y se documentan las diferencias.
3. **PLAN.md Fase 1 "función única `get_main_stylesheet(theme)`"**: en el
   código real existen dos generadores y dos constantes consumidos por
   `main_window.py` y `toolbar.py`. La fachada conserva ambos generadores y
   las constantes (compatibilidad), no una función única.
4. **El tema no se persistía hoy**: `shell/assembly.py` hardcodea
   `current_theme = "light"`; `PreferencesDialog` ofrece claro/oscuro pero la
   elección se perdía al cerrar. La persistencia `ui/theme` (M21) se añade en
   esta fase como parte de la fundación.
5. **`icons.py` sin colores propios es objetivo de Fase 2**, no de esta fase
   (instrucción explícita: no migrar iconos aún; `icons.py` intacto).
6. **Entorno de tests**: el worktree no tiene `.venv`; el venv del checkout
   principal (`/home/unison-pjgv/Documentos/GitHub/Chemuson/.venv`, Python
   3.11.16) trae el stack runtime (PyQt6 6.11, rdkit, numpy, Pillow) pero no
   pytest/ruff. Se usa entorno efímero `uv run --no-project --offline` sobre
   el intérprete del venv con pytest 9.1.1/ruff resueltos desde la caché de
   uv (sin instalaciones nuevas ni dependencias adicionales).
