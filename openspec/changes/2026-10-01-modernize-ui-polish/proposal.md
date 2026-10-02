# Proposal: Pulido final de la UI moderna (Fase 7)

## Why

Las fases 1–6 están aprobadas manualmente en KDE/Wayland. Quedan pendientes
los puntos de "polish" de `docs/ui-modernization/PLAN.md` §Fase 7 y residuos
concretos detectados por auditoría (ver `design.md`):

- **Contraste**: el token light `text3` (`#94A3B8`) tiene ratio ~2.3:1 sobre
  `#F1F5F9` (bajo WCAG AA); es el token de texto terciario/deshabilitado.
- **Disabled states**: `:disabled` en QSS solo cambia color a `text3` sin
  opacidad ni fondo; un botón deshabilitado (p. ej. Deshacer) no se distingue
  inequívocamente de un botón en reposo.
- **Tooltips**: convención esperada `Nombre` / `Nombre (Shortcut)` tomada de
  `QAction.shortcut()`; hoy el rail usa strings hardcodeados en `_RailSpec` y
  la pestaña "Nuevo" declara `(Ctrl+N)` sin referencia a la `QAction`.
- **HiDPI**: el `IconProvider` ya es DPR-aware, pero el thumbnail de
  plantillas rasteriza `QPixmap(88, 56)` fijo (píxeles físicos) con pen/font
  hardcodeados: a 125/150/200 % queda difuminado.
- **Onboarding**: no existe overlay de primera ejecución (solo el diálogo
  legado `QuickStartDialog`). El plan pide 3 puntos (rail, lienzo, panel
  derecho) con "No mostrar de nuevo" persistido.
- **Iconos provisionales**: la barra de formato de texto usa glifos
  tipográficos (`≡`, `☰`, …) para alineación/sub/superíndice.
- **Referencias stale**: docstrings/comentarios de `command_palette.py`,
  `qss.py`, `tokens.py`, `main_window.py` y `app_bar.py` aún dicen "Ctrl+K"
  para la paleta (hoy `Ctrl+P`; `Ctrl+K` es Limpiar 2D 1 paso).
- **Paleta/Plantillas**: el registro se reconstruye al abrir la paleta, pero
  si la paleta está abierta y cambia la biblioteca, el registro queda stale.
- **Manual**: `docs/MANUAL_USUARIO.md` §2/§3/§6.5/§9/§10/§11/§12 describen la
  UI antigua (barra principal, barras izquierda/derecha) y no documentan la
  paleta ni los atajos reales.

## What Changes

- `theme/tokens.py`: subir `text3` light a `#5E6E82` (AA en `#F1F5F9`/`#FFFFFF`/`#F8FAFC`); dark intacto.
- `theme/qss.py`: `:disabled` inequívoco para `QToolButton`, `#railBtn`,
  celdas `flyItem`/`paletteItem` (opacidad + fondo + color de etiqueta);
  orden de reglas garantiza que `:disabled` gana sobre `:hover`.
- `tool_rail.py`: tooltips construidos desde `QAction.shortcut()` cuando
  existe; sin duplicar la infraestructura de presentación.
- `document_tabs.py`: tooltip de "Nuevo" referenciado a `action_new`.
- `template_browser_service.py`: thumbnail DPR-aware (márgenes/pen/font
  escalados) + marco uniforme (solo presentación; sin tocar estructura
  química del grafo).
- `text_toolbar.py`: 5 SVG nuevos (`i-align-left/center/right`,
  `i-subscript`, `i-superscript`) sustituyen glifos provisionales; `B/I/U`
  se conservan (convención de editores).
- Nuevo `gui/onboarding.py`: `OnboardingOverlay` (3 pasos, Anterior/Siguiente/
  Cerrar, "No mostrar de nuevo") inyectado en `shell/assembly.py` solo en
  primera ejecución; persistencia `platform.settings` clave
  `ui/onboarding/completed`. No toca escena/canvas internals.
- `main_window.py`: `_refresh_template_views` refresca también el registro de
  la paleta (arreglo pequeño, sin tocar `TemplateLibrary` ni química).
- Comentarios/docstrings stale: `Ctrl+K` → `Ctrl+P` donde la paleta es el
  sujeto (`command_palette.py`, `qss.py`, `tokens.py`, `main_window.py`,
  `app_bar.py`).
- `docs/MANUAL_USUARIO.md`: reescribir §2/§3/§6.5/§10/§11/§12 para la UI
  moderna (AppBar, DocumentTabs, ToolRail + flyouts, SidePanel,
  CommandPalette, tema claro/oscuro, onboarding) y atajos reales
  (`Ctrl+P` paleta, `Ctrl+K` Clean2D 1 paso, `Ctrl+Shift+K` publicación,
  `Ctrl+Alt+K` conformero).

## Impact

- Archivos: `theme/tokens.py`, `theme/qss.py`, `tool_rail.py`,
  `document_tabs.py`, `template_browser_service.py`, `text_toolbar.py`,
  `theme/icons/` (+5 SVG), nuevo `onboarding.py`, `shell/assembly.py`,
  `main_window.py`, `command_palette.py`, `app_bar.py`,
  `docs/MANUAL_USUARIO.md`, tests nuevos.
- Sin cambios en `clean2d/`, `chemname/`, `chemio/persistence.py`; sin tocar
  `gui/canvas/`, `gui/editor2d/`, geometría, selección, `tool_id`, undo
  semantics ni química de plantillas.
- Tests: suite verde idéntica al baseline (1813 passed / 20 skipped /
  5 failed conocidos) + tests nuevos de tooltips, onboarding, HiDPI,
  disabled, DPR de thumbnails y refresh del registro.
