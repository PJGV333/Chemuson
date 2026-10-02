# Tasks: Pulido final de la UI moderna (Fase 7)

## 1. OpenSpec y baseline

- [x] 1.1 Confirmar rama `ui/modernization`, HEAD local/remoto `50a8241`
  ("Restore Clean2D shortcut and finalize command palette") y árbol limpio.
- [x] 1.2 Auditoría amplia de polish (delegada a Ornith 9B vía mecanismo de
  delegación de Hermes; verificado punto por punto por Qwen contra el código):
  contraste light `text3`, `:disabled` débil, tooltips ad-hoc, thumbnails
  fijos 88×56, glifos provisionales, comentarios stale Ctrl+K, onboarding
  ausente, staleness de plantillas en paleta, manual desactualizado.
- [x] 1.3 Ejecutar compileall, colección de pytest, suite completa y Ruff
  scoped; registrar en `baseline.md` (logs en `/tmp/f7_*.log`).
- [x] 1.4 Validar este cambio con
  `openspec validate 2026-10-01-modernize-ui-polish --strict` **antes** de
  tocar implementación.

## 2. Contraste y disabled (tokens.py, qss.py)

- [x] 2.1 `text3` light → `#5E6E82` (AA en `#F1F5F9`/`#FFFFFF`/`#F8FAFC`); bloque dark intacto.
- [x] 2.2 `:disabled` inequívoco: `QToolButton`, `#railBtn`,
  `flyItem` (+`#flyLbl`), `paletteItem`, `QPushButton`/flat — opacidad +
  fondo/color; `:disabled` tras `:hover` en la hoja (gana).
- [x] 2.3 Sin lógica de enabled en widgets: la `QAction` sigue siendo la
  fuente de verdad (undo/redo intactos).

## 3. Tooltips uniformes

- [x] 3.1 Rail: tooltip desde `QAction.shortcut().toString()` cuando existe;
  en caso contrario, tooltip del spec (letras contextuales V/B/L/R/C/T/N/G/E/O).
- [x] 3.2 `document_tabs.py`: tooltip de "Nuevo" referenciado a
  `action_new.shortcut()`.
- [x] 3.3 `SearchPill`: tooltip referenciado a la `QAction` de apertura
  ("Buscar o ejecutar… (Ctrl+P)").
- [x] 3.4 Sin duplicar infraestructura de presentación (reutilizar patrón de
  `command_palette._action_tooltip`): helpers locales `_tooltip_from_action`/`_tooltip_for_action` siguen la misma convención `Nombre (Shortcut)`. Nota: el tooltip se fija en la `QAction` para sobrevivir a los re-sync de Qt.

## 4. HiDPI y thumbnails

- [x] 4.1 `template_preview_icon`: pixmap `88×56` × dpr con
  `setDevicePixelRatio(dpr)`; margen/pen/font escalados; grafo intacto.
- [ ] 4.2 Smoke offscreen `QT_SCALE_FACTOR=2`: ventana, rail, app bar,
  side panel; iconos no nulos a dpr 1 y 2.

## 5. Onboarding (nuevo gui/onboarding.py)

- [x] 5.1 `OnboardingOverlay`: 3 pasos (Tool Rail / Canvas / SidePanel) con
  textos exactos del brief; Anterior/Siguiente/Cerrar + "No mostrar de nuevo".
- [x] 5.2 Inyección en `shell/assembly.py` solo en primera ejecución vía
  `platform.settings` clave `ui/onboarding/completed`; sin tocar escena/canvas.
- [x] 5.3 Tests: 3 pasos, persistencia, no se repite, "no mostrar" persiste.

## 6. Iconos provisionales (text_toolbar.py)

- [x] 6.1 SVGs `i-align-left/center/justify`, `i-subscript`, `i-superscript`
  (24 px, trazo 1.75) + integración vía `IconProvider` (tint por tema, DPR).
- [x] 6.2 Botones de alineación/sub/sup usan los SVG; `B/I/U` se conservan.

## 7. Plantillas en CommandPalette (arreglo pequeño)

- [x] 7.1 `_refresh_template_views` refresca también el registro de la paleta
  (`_refresh_command_registry_templates`); sin tocar `TemplateLibrary`/química.
- [x] 7.2 Test: plantilla creada/eliminada → registro coherente al abrir.

## 8. Referencias stale

- [x] 8.1 `Ctrl+K` → `Ctrl+P` donde la paleta es el sujeto:
  `command_palette.py:1`, `qss.py:970`, `tokens.py:167`,
  `main_window.py:312/316-317`, `app_bar.py:131`.
- [x] 8.2 Conservar `structure_actions.py`, `assembly.py`, `tool_rail.py:396`
  (ahí `Ctrl+K` es correctamente Clean2D) y PLAN.md/mockup (históricos).

## 9. Manual (docs/MANUAL_USUARIO.md)

- [x] 9.1 §2: mapa moderno (AppBar, DocumentTabs, ToolRail+flyouts,
  SidePanel, texto de formato, lienzo, estado, tema claro/oscuro, onboarding).
- [x] 9.2 §3: referencias al panel de herramientas (no "barras izq/der").
- [x] 9.3 §6.5: paneles como tabs del SidePanel (menú Ver back-compat).
- [x] 9.4 §10: "Barra de aplicación" + subsección "Paleta de comandos
  (Ctrl+P)" con atajos reales (Ctrl+K, Ctrl+Shift+K, Ctrl+Alt+K).
- [x] 9.5 §11/§12: renombrar a ToolRail y a herramientas de anotación/símbolos
  (flyouts del mismo rail); contenido conservado; §9 gana la ruta rápida de la
  paleta. Sin inventar atajos.

## 10. Evidencia y validación final

- [x] 10.1 Screenshots reales (offscreen): LIGHT y DARK × (1440×900, 980×600)
  + onboarding + Templates + CommandPalette; HiDPI 200 % (y 125/150 % si
  estable).
- [x] 10.2 Revisión final del diff con Ornith (scope creep, imports,
  duplicación, shortcut regressions, theme inconsistencies, regresiones
  visuales; máx. 20 hallazgos); Qwen decide cada hallazgo. Veredicto
  PROCEDE CON OBSERVACIONES; resueltos: propuesta `#5E6E82` (no `#64748B`),
  spec/design `i-align-justify` (no `i-align-right`), etiqueta onboarding
  "Siguiente" en todos los pasos (spec estricto), y limitación HiDPI
  `QT_SCALE_FACTOR` documentada en design.md.
- [x] 10.3 Tests nuevos + suite completa + compileall + Ruff scoped +
  `git diff --check` + OpenSpec strict; comparar contra baseline de esta PC.
  Suite: 5 failed (idénticos a los 5 preexistentes), 1830 passed, 20 skipped
  (+17 nuevos vs. baseline 1813); Ruff: solo el error preexistente
  (F401 `math`); `git diff --check` OK; OpenSpec strict válido.
- [x] 10.4 Commit `Polish modern UI and add onboarding` (+
  `Update modern UI user guide` si el manual va separado); push sin
  `--force`; verificar HEAD local/remoto y worktree limpio. Hecho:
  commit 22496d9 (polish) + 1238263 (manual); push origin/ui/modernization;
  local == remoto == 1238263b8d14657e29c5d1166c3b7bf3d8270ec5; worktree limpio.
