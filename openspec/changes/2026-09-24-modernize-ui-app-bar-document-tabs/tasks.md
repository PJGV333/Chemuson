# Tasks: Fase 3 — App bar y pestañas de documento

## OpenSpec
- [x] 1. Capturar baseline (git status, compileall, collect, pytest, ruff) en `baseline.md`.
- [x] 2. `openspec validate 2026-09-24-modernize-ui-app-bar-document-tabs --strict` ANTES de implementar.

## Iconos (Fase 2, extensión)
- [x] 3. Añadir 8 SVG nuevos en `src/chemuson/gui/theme/icons/` (plus, search, moon, sun, sliders, flask, x, doc) con el contrato 24×24/`currentColor`; actualizar `LICENSE.txt` si procede.

## Componentes nuevos
- [x] 4. `src/chemuson/gui/document_tabs.py`: `DocumentTabBar` (QTabBar espejo: elide right, scroll, movable, icono doc, punto de suciedad, botón cerrar, corner `+`, señales `tabActivated`/`newDocumentRequested`, `sync_tabs`/`set_tab`/`select`, `refresh_icons(theme_name)`).
- [x] 5. `src/chemuson/gui/app_bar.py`: `SearchPill` (placeholder sin atajo) y `AppBar` (54 px; marca flask+nombre+versión; undo/redo/tema/ajustes con `setDefaultAction`; `refresh_icons(theme_name)`; `sync_tabs`/`select`/`set_tab`).
- [x] 6. `theme/qss.py`: sección de QSS de tokens para la app bar y las pestañas.
- [x] 7. `theme/tokens.py`: métrica `appbarH: 54` en `METRICS`.

## Integración
- [x] 8. `tab_manager.py`: kwarg opcional `on_change` en `CanvasTabManager` (invoke en `create_document_tab`/`discard_canvas`); comportamiento por defecto idéntico.
- [x] 9. `main_window_ui_builder.py`: crear `action_preferences` y `action_theme_toggle` en `create_local_actions`; `build_menu_bar` deja de crear `action_preferences` (mismo objeto).
- [x] 10. `shell/assembly.py`: wrapper central `[app_bar, tabs]`, `tabs.tabBar().hide()`, observer + conexiones (tabActivated→setCurrentIndex, closeRequested→`_on_tab_close_requested`, newDocumentRequested→`action_new.trigger()`, tabMoved→`moveTab`, callback `on_tab_updated` (de `update_tab_title`)→resync de texto/suciedad, `currentChanged`→select), `main_toolbar` oculta.
- [x] 11. `main_window.py`: `_sync_app_bar_tabs` (resync completo + por pestaña), `select` desde `_on_tab_changed`, `action_theme_toggle` sync en `_apply_theme`, `app_bar.refresh_icons(resolved)` en `_apply_theme`.
- [x] 12. `architecture/modules.yml`: registrar nuevos módulos/tests en M08.

## Verificación
- [x] 13. `tests/test_ui_app_bar_tabs.py`: app bar/tabs unitarios + ventana real (inicial, nuevo, suciedad, cierre, cambio, reorder, QAction compartidas, enabled/disabled, tema light→dark→light, resize 1440×900/980×600, sin atajo en pill).
- [x] 14. Extender `tests/test_ui_svg_icons.py`: 8 SVG nuevos en el inventario (63) + resolución del provider.
- [x] 15. Smoke Qt offscreen de la ventana real (light/dark, varias pestañas, dirty, undo enabled/disabled).
- [x] 16. Suite completa + ruff scoped + `git diff --check` + `compileall` (sin regresiones vs baseline).
- [x] 17. Capturas `docs/ui-modernization/app-bar-phase-shots/` (light/dark 1440×900 y 980×600, varias pestañas, pestaña modificada, undo/redo estados) + script reproducible + README de comparación.
- [x] 18. Revalidar OpenSpec `--strict` si cambió algún doc; marcar tasks completas.
