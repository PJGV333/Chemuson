# Tasks: Paleta de comandos Ctrl+K (Fase 6)

## 1. OpenSpec y baseline

- [x] 1.1 Confirmar rama `ui/modernization`, HEAD local/remoto `02c8ea8` y árbol limpio.
- [x] 1.2 Revisar PLAN.md §2.4/Fase 6, mockup-ui.html, spike `palette.py`, app_bar, ui_builder, actions/, shell, side_panel, template_browser_service, settings, modules.yml y tests de AppBar/shortcuts/docks.
- [x] 1.3 Ejecutar compileall, colección de pytest, suite completa y Ruff; registrar en `baseline.md`.
- [x] 1.4 Validar este cambio con `openspec validate 2026-10-01-modernize-ui-command-palette --strict` **antes** de tocar implementación.

## 2. CommandPalette (command_palette.py)

- [x] 2.1 Implementar `CommandEntry` (action, section, keywords, icon) y `CommandRegistry` con deduplicación por `id(action)` y `register(action, section, keywords)`.
- [x] 2.2 Implementar `CommandPalette` (overlay centrado, input, lista filtrada) con ranking prefijo>substring determinista.
- [x] 2.3 Implementar teclado: `↑`/`↓`, Enter (1 vez + cierra), Esc (cierra sin ejecutar), clic ejecuta; fila disabled no ejecutable; checkable por delegación a `triggered`.
- [x] 2.4 "Último comando usado" local a sesión (sin persistir) como fila 0 con query vacío.
- [x] 2.5 Supresión de atajos de herramienta mientras el foco está en el `QLineEdit` (delegación al dispatcher existente).

## 3. Fuentes de comandos y montaje

- [x] 3.1 Registrar las `QAction` reales de Archivo, Editar, Ver, Estructura, Análisis, exportaciones (PNG/SVG/PDF/CML/SMILES), canvas/view, preferencias/tema.
- [x] 3.2 Reutilizar `window.side_panel_actions` (7 páginas) sin llamar `show_page()` en paralelo.
- [x] 3.3 Plantillas: reutilizar las `QAction` del menú dinámico `window.templates_menu` (contrato existente) vía `refresh_templates()`; adaptador mínimo documentado como fallback (design D5). No tocar contenido químico.
- [x] 3.4 Construir una única `CommandPalette` en `shell/assembly.py` y una única `QAction` de apertura `action_command_palette` (Ctrl+K, WindowShortcut, `window.addAction`).
- [x] 3.5 Conectar `app_bar.search_pill.activated` → `action_command_palette.trigger` (mismo camino que Ctrl+K).

## 4. Migración de Ctrl+K

- [x] 4.1 Retirar `Ctrl+K` de `action_clean_2d_full` (quitar shortcut + `window.addAction`), conservando `QAction`, handler y función.
- [x] 4.2 No asignar otro shortcut a `action_clean_2d_full` en esta fase; `Ctrl+Shift+K`/`Ctrl+Alt+K` intactos.
- [x] 4.3 Actualizar únicamente tests/documentación que afirmaban `Ctrl+K → Clean2D quick`.

## 5. SearchPill del AppBar (mínimo)

- [x] 5.1 Mostrar el badge `Ctrl K` y quitar "(próximamente)" del tooltip.
- [x] 5.2 Convertir la `SearchPill` en entrada real (clic → misma `QAction` de apertura) sin convertirla en editor permanente ni rediseñar el AppBar.

## 6. Tokens/QSS y arquitectura

- [x] 6.1 Añadir métrica `paletteW` (560) y QSS tokenizado del overlay/tarjeta/input/sections/rows/selected (light+dark) sin colores fuera de tokens.
- [x] 6.2 Registrar `command_palette.py` en M08 de `architecture/modules.yml`.
- [x] 6.3 Añadir contrato AST: `command_palette.py` NO importa `clean2d`/`chemname`/`chemio.persistence`/`gui.canvas`/controllers químicos.

## 7. Tests obligatorios

- [x] 7.1 `tests/test_command_palette.py`: ventana construye una única paleta; registro ≥60 únicos; sin duplicar `QAction`.
- [x] 7.2 Las 7 páginas del SidePanel buscables; PNG/SVG/PDF/CML/SMILES buscables; Clean2D quick buscable por su `QAction` histórica.
- [x] 7.3 Ctrl+K no pertenece a `action_clean_2d_full`; Ctrl+K abre la paleta; clic en SearchPill abre la misma ruta; hint `Ctrl K` visible.
- [x] 7.4 Ctrl+Shift+K y Ctrl+Alt+K conservan su comportamiento.
- [x] 7.5 Prefijo antes que substring; matching por keywords; ↑/↓ cambia selección; Enter dispara 1 vez; Esc cierra sin disparar; disabled no ejecuta; checkable conserva semántica; al ejecutar se cierra; reabrir no duplica; light→dark→light; 980×600 dentro de ventana; tool shortcuts no interfieren al escribir; acciones siguen accesibles por menú; sin conflicto de Ctrl+K.
- [x] 7.6 Comandos de diálogos modales: probar wiring/trigger con `QAction` segura o mocking puntual (sin bloquear offscreen).

## 8. Evidencia y validación final

- [x] 8.1 Capturas offscreen con la ventana real: light/dark × 1440×900 (query vacío, "export", "valid") y 980×600.
- [ ] 8.2 Tests dirigidos + AppBar + shortcuts + SidePanel + architecture + compileall + Ruff scoped + `git diff --check` + OpenSpec strict + suite completa; comparar contra baseline de esta PC.
- [ ] 8.3 Revisar diff completo; commit `Add command palette`; push sin `--force`; verificar HEAD local/remoto y worktree limpio.
