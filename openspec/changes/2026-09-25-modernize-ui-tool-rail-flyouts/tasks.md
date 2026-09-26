# Tasks: Fase 4 — Rail de herramientas unificado + flyouts

## OpenSpec
- [x] 1. Capturar baseline (git status, compileall, collect, pytest, ruff + inventario 1:1 de tool_ids/acciones/señales/atajos) en `baseline.md`.
- [x] 2. `openspec validate 2026-09-25-modernize-ui-tool-rail-flyouts --strict` ANTES de implementar.

## Componentes nuevos
- [x] 3. `src/chemuson/gui/flyout.py`: `Flyout` + `FlyoutItem` + `FlyoutFooter` + `FlyoutCell` (244 px, cabecera título+Kbd Esc, grid N columnas, pie ≤3 botones, `show_near`, cierre Esc/clic-fuera/selección, señal `closed`, `on_select`, `checked`/`active` por celda).
- [x] 4. `src/chemuson/gui/tool_rail.py`: `ToolRail` (58 px) + `ToolRailButton` (icono 21 px + kbd-hint + activo) + `ToolShortcutDispatcher` (event filter contextual). Construido desde `ChemusonToolbar` + `SymbolPaletteToolbar` (D1: delegación 1:1; introspección de `QMenu` originales; flyouts por grupo D6; `set_active_tool`/`clear_active`/`refresh_icons`).
- [x] 5. `theme/tokens.py`: `METRICS` += `railW: 58`, `flyoutW: 244`.
- [x] 6. `theme/qss.py`: QSS de tokens para `#toolRail`, `#railBtn`, `#railKbd`, `#railSep`, `#flyout`, `#flyLbl`, `#flyoutTitle`, `#flyFoot`, `#kbdPill`.

## Integración
- [x] 7. `shell/assembly.py`: montar `ToolRail` (QToolBar wrapper en LeftToolBarArea), ocultar `toolbar` y `symbols_toolbar` (sin eliminar), conectar `tool_changed` → `tool_rail.set_active_tool`, dispatcher de atajos (V/A/L/B/R/C/T/N/G/E/O), `_apply_theme` refresca el rail tras los toolbars, `_clear_active_tool_selection` → `tool_rail.clear_active()`.
- [x] 8. `architecture/modules.yml`: M08 paths += `tool_rail.py`, `flyout.py`; internal_api += `tool_rail`, `flyout`; nota de fase 4.

## Verificación
- [x] 9. `tests/test_ui_tool_rail.py`: inventario 1:1 (contadores de celdas vs menús originales), delegación (flyout → señales del toolbar), footers (ring size, periodic, diálogos, presets), estado activo (tool_changed, clear en cambio de pestaña, exclusión del `QActionGroup`, sin doble emisión), atajos contextuales (foco canvas / QLineEdit / modificadores / diálogo modal / Ctrl+K intacto), toolbars ocultos y presentes, `QMenuBar` visible, tema light→dark→light, métricas.
- [x] 10. Smoke Qt offscreen de la ventana real (rail visible, flyout abierto/cerrado, atajos, sin tracebacks).
- [x] 11. Suite completa + ruff scoped + `git diff --check` + `compileall` (sin regresiones vs baseline; los 4 fallos conocidos se mantienen sin tocar).
- [x] 12. Capturas `docs/ui-modernization/tool-rail-phase-shots/` (light/dark 1440×900, rail, flyouts de enlace/símbolos/orbitales) + script reproducible + README.
- [x] 13. Revalidar OpenSpec `--strict` si cambió algún doc; marcar tasks completas.
- [x] 14. Commit `Add unified tool rail and flyouts` + push a `origin ui/modernization` (si la autenticación lo permite; si no, dejar el commit local y reportar).
