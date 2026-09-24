# Design: UI Theme Foundation and Design Tokens

## Resumen

`src/chemuson/gui/theme/` se convierte en la única fuente de verdad visual de
la UI: tokens (light/dark) → generadores QSS + QPalette + fuente. `styles.py`
queda como fachada de compatibilidad. El tema se aplica desde un punto central
(`theme.apply_theme`) que la ventana actual ya invoca al final del ensamblaje
y en cada cambio de tema.

## Archivos y módulos afectados

| Archivo | Cambio | Módulo |
|---|---|---|
| `src/chemuson/gui/theme/__init__.py` | Nuevo: API pública (`apply_theme`, `resolve_theme_name`, `set_theme_from_system`, re-exports) | M08 |
| `src/chemuson/gui/theme/tokens.py` | Nuevo: tabla light/dark + métricas (spacing/radios/tipografía) | M08 |
| `src/chemuson/gui/theme/qss.py` | Nuevo: `get_main_stylesheet` / `get_tool_palette_stylesheet` desde tokens | M08 |
| `src/chemuson/gui/theme/palette.py` | Nuevo: `build_qpalette`, `theme_font` | M08 |
| `src/chemuson/gui/theme/icon_provider.py` | Nuevo: `IconProvider` SVG→QIcon/QPixmap, caché + HiDPI (Fase 2 prepara) | M08 |
| `src/chemuson/gui/styles.py` | Reescrito como fachada: re-exports + alias legados | M08 |
| `src/chemuson/gui/main_window.py` | `_apply_theme` delega en `theme.apply_theme`; carga/persistencia de `ui/theme` | M08 |
| `src/chemuson/gui/shell/assembly.py` | Carga del tema persistido tras crear `self._settings` | M08 |
| `src/chemuson/platform/settings.py` | `UiPreferences`, `load_ui_preferences`, `save_ui_preferences` | M21 |
| `src/chemuson/platform/__init__.py` | Export de la nueva API UI (aditivo) | M21 |
| `architecture/modules.yml` | M08: paths/`internal_api`/tests de `theme/`; M21: `public_api` nuevo | — |
| `tests/test_ui_theme_foundation.py` | Nuevo: tests de la fundación | M08 |
| `tests/test_platform_settings.py` | Extendido: round-trip `ui/theme` | M21 |

**No se tocan**: `gui/icons.py` (intacto, con `set_icon_theme`), `gui/canvas/`
(M09), `clean2d/`, `chemname/`, `chemio/persistence.py`, `gui/style.py`
(`DrawingStyle` es parámetro de dibujo químico, no tema UI), `gui/toolbar.py`,
`gui/docks.py`, `gui/text_toolbar.py`, `gui/dialogs/` (el combo de tema ya
funciona y sigue haciéndolo; solo se persiste la elección).

## Decisiones de diseño

### D1. Tokens: el spike es la referencia

Los valores y el vocabulario de tokens se copian de
`docs/ui-modernization/pyqt6-spike/theme.py` (aprobado visualmente, commit
`59e977d`), no de la tabla textual de PLAN.md §2.2 (ver discrepancias 1–2 del
propuesto). Se añade la capa de métricas que el spike no centralizaba:
spacing (grilla 4/8/12/16 px), radios (10/9/8) y tipografía (13/12/11 px)
según PLAN.md §2.2. Los tokens químicos (CPK, hoja blanca del canvas,
`sheet`/`sheetGrid`/`canvasBg`) se incluyen pero **no se aplican** al canvas
en esta fase (el canvas conserva su render propio; M09 intacto).

### D2. Un solo generador de QSS por hoja

`qss.py` expone exactamente los dos contratos existentes
(`get_main_stylesheet`, `get_tool_palette_stylesheet`) construyendo las
reglas por f-string sobre `get_tokens(theme)`. La estructura de reglas se
migra de `styles.py` actual (menubar, menús, toolbars, toolbuttons, docks,
statusbar, scrollbars, tablas, árboles, labels, diálogos, botones,
lineedits, combos, spins, checkboxes, radios, tabs, groupboxes, tooltips,
`#palette_grid`) con valores del spike: superficies claras/limpias (`bg`,
`surface`, `surface2`), cyan de acento (`accent`), bordes suaves (`border`),
botones modernos y estados hover/pressed/checked/disabled coherentes
(hover `surface2`, pressed `surface3`, checked `accentSoft`+`accentBorder`,
disabled `text3`). No se usan propiedades CSS inexistentes en Qt (sin
`box-shadow` ni `transitions`; `opacity` QSS no se usa: los estados
disabled se pintan con colores explícitos de tokens).

### D3. QPalette además de QSS

`build_qpalette(theme)` mapea tokens a roles (Window=`bg`, Base=`surface2`,
AlternateBase=`surface3`, Button=`surface`, Text/WindowText/ButtonText=`text1`,
PlaceholderText=`text3`, Highlight=`accent`, HighlightedText=`onAccent`,
ToolTipBase/ToolTipText=`surface`/`text1`). Hoy no se aplica QPalette (solo
QSS); aplicarlo alinea widgets que heredan del sistema (menús nativos,
popups) con los tokens, igual que hace el spike.

### D4. `styles.py` como fachada

`chemuson.gui.styles` re-exporta desde `chemuson.gui.theme`:
`get_main_stylesheet`, `get_tool_palette_stylesheet`, `DEFAULT_THEME`
(=`light`), `MAIN_STYLESHEET`, `TOOL_PALETTE_STYLESHEET` (generadas en
import con el tema por defecto, como hoy) y `LIGHT_COLORS`/`DARK_COLORS`
como alias de tokens (las 22 claves legadas mapeadas 1:1 a tokens; se marca
`DeprecationWarning` en el docstring, no se emite warning en runtime para no
ensuciar consolas). Callers reales: `main_window.py` (los dos generadores) y
`toolbar.py` (`TOOL_PALETTE_STYLESHEET`); no hay uso de las paletas legadas.

### D5. Punto central de aplicación

`theme.apply_theme(target, theme_name)`:
1. `theme_name_resuelto = resolve_theme_name(theme_name)` (normaliza y
   resuelve `"system"` → light/dark con `QStyleHints.colorScheme()` si
   existe; cualquier valor desconocido → `"light"`, sin excepciones).
2. `target.setFont(theme_font())`
3. `target.setPalette(build_qpalette(resuelto))`
4. `target.setStyleSheet(get_main_stylesheet(resuelto))`

Acepta `QApplication` o `QWidget` (la ventana actual). `ChemusonWindow._apply_theme`
llama a esto y conserva después el QSS de paletas para `toolbar`/
`symbols_toolbar` y el refresco de iconos (comportamiento actual).
`set_theme_from_system(target)` = `apply_theme(target, "system")`
(una-shot; no hay señal portable de cambio de esquema en Qt → el
re-apply llega con la Fase de ajustes si procede).

### D6. Persistencia mínima en M21

`UiPreferences(theme: str = "light")` + load/save sobre la clave `ui/theme`
(valores admitidos: `light`, `dark`, `system`; cualquier otro → `light`).
M21 no importa QtWidgets (solo `QSettings`/Protocol, como hoy).
`assembly.py` carga tras crear `self._settings` (el `current_theme = "light"`
inicial sigue como fallback antes de que exista `_settings`);
`_apply_preferences` persiste al cambiar. `toggle_theme` y el combo de
Preferencias siguen operando solo light/dark (el valor `system` es
aceptado/normalizado por la API, no se expone en la UI esta fase).

### D7. IconProvider: infraestructura, no migración

`icon_provider.IconProvider` (modelo del spike, adaptado al paquete):
- `icon(name, color, size=20) -> QIcon` y `pixmap(name, color, size=20) -> QPixmap`
- tinte por sustitución de `currentColor` en el SVG (sin iterar píxeles)
- caché de bytes por nombre y caché de `QIcon`/`QPixmap` por
  `(name, color, size)` (theme-aware: el color del token es la clave)
- HiDPI: `devicePixelRatio` configurable (dpr del `QGuiApplication`)
- SVG estáticos desde `src/chemuson/gui/theme/icons/` (carpeta que la Fase 2
  poblará; un nombre ausente devuelve `QIcon()`/`QPixmap()` vacíos: fallo
  visible, sin excepción) y glifos dinámicos `glyph:<label>|<shape>|<size>`
  (mismo esquema del spike) para probar el pipeline sin assets.
- `icons.py` NO se modifica; la delegación `icons.py → IconProvider` es Fase 2.

### D8. Registro arquitectónico

M08: `paths` += `src/chemuson/gui/theme/`; `internal_api` += `theme`;
`tests` += `tests/test_ui_theme_foundation.py`. M21: `public_api` +=
`UiPreferences`, `load_ui_preferences`, `save_ui_preferences`. Sin
dependencias nuevas (el subpaquete solo importa PyQt6 ya dependiente de M08;
M21 no gana imports). Sin excepciones temporales ni circularidades.

## Riesgos de regresión y mitigación

| Riesgo | Mitigación |
|---|---|
| QPalette nuevo cambia widgets que heredaban del sistema | Baseline de suite completa antes/después; smoke visual offscreen light/dark de la ventana real; las reglas QSS cubren los widgets principales |
| Fuentes: ancho de widgets calculado con la fuente anterior (pitfall del spike) | El tema (fuente incluida) se aplica al final del ensamblaje, igual que hoy; los widgets se crean antes que el `_apply_theme` final, idéntico al orden actual |
| Rotura de callers de `styles.py` | Fachada 1:1 (mismos nombres); test explícito de compatibilidad; grep verificó los únicos 2 call sites |
| QSS con propiedades no soportadas por Qt | Solo propiedades QSS nativas (ver D2); validación visual con smoke offscreen |
| Canvas (M09) alterado | `git diff` restringido a los archivos del design; la suite de canvas/Clean2D de la baseline debe pasar igual |
| M21 gana acoplamiento GUI | `UiPreferences` es un dataclass de string sin imports de Qt-GUI; test de round-trip con `FakeSettings` (sin QApplication) |

## Pruebas

- `tests/test_ui_theme_foundation.py`: tokens light/dark resueltos y
  diferenciales; todos los colores de tokens válidos (`QColor.isValid`);
  métricas presentes; `resolve_theme_name` (incl. `system` y valores basura);
  `apply_theme` light/dark sin excepciones sobre QApplication y QMainWindow
  (offscreen); QSS contiene valores de tokens del tema; QPalette coherente
  con tokens; compatibilidad de la fachada `styles.py` (nombres + constantes);
  ventana real `ChemusonWindow` sin excepciones con ambos temas
  (`toggle_theme`); IconProvider: caché (misma instancia por clave),
  pixmap no nulo con tamaño/dpr correctos, glifos dinámicos, nombre ausente
  → vacío sin excepción.
- `tests/test_platform_settings.py` (extendido): load normaliza valores
  inválidos a `light`, conserva `dark`/`system`, save escribe `ui/theme`.
- Smoke Qt offscreen: creación de `ChemusonWindow` + aplicación de ambos
  temas + capturas de referencia en
  `docs/ui-modernization/foundation-shots/`.
- Baseline/verificación (AGENTS.md §1.2/§1.3): `git status --short`,
  `compileall`, `pytest --collect-only -q`, `pytest -q`, `ruff check ...
  --select F401,F811,F821,E722,E741` registrados en `baseline.md` antes y
  después.
