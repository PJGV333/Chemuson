# Tasks: UI Theme Foundation and Design Tokens

## 1. Preparación y baseline (Fase 0)

- [x] Verificar rama `ui/modernization`, árbol limpio y `git pull --ff-only`.
- [x] Auditoría de arquitectura: ubicar estilos/temas, iconos, settings,
      shell y tests de GUI; registrar discrepancias vs PLAN.md (ver
      proposal.md §Discrepancias).
- [x] Localizar entorno de test correcto (uv efímero offline sobre el venv
      del checkout principal; sin instalaciones nuevas).
- [x] Capturar baseline completa en `baseline.md` (git status, compileall,
      pytest collect + suite completa, ruff scoped).
- [x] Validar el cambio OpenSpec en modo estricto antes de implementar.

## 2. Design tokens (Fase 1)

- [x] `src/chemuson/gui/theme/tokens.py`: tabla light/dark del spike aprobado
      + métricas (spacing/radios/tipografía) + `get_tokens`,
      `theme_color`, `THEME_NAMES`, `DEFAULT_THEME_NAME`.
- [x] `src/chemuson/gui/theme/palette.py`: `build_qpalette(theme)` y
      `theme_font(...)` desde tokens.

## 3. QSS + punto central de aplicación (Fase 1)

- [x] `src/chemuson/gui/theme/qss.py`: `get_main_stylesheet` /
      `get_tool_palette_stylesheet` 100 % desde tokens (reglas actuales
      migradas, look del spike: superficies limpias, cyan accent, bordes
      suaves, hover/pressed/checked/disabled coherentes).
- [x] `src/chemuson/gui/theme/__init__.py`: `apply_theme`,
      `resolve_theme_name` (con `system` vía QStyleHints),
      `set_theme_from_system`, re-exports.
- [x] `src/chemuson/gui/styles.py`: fachada de compatibilidad 1:1
      (generadores, `DEFAULT_THEME`, `MAIN_STYLESHEET`,
      `TOOL_PALETTE_STYLESHEET`, `LIGHT_COLORS`/`DARK_COLORS` como alias de
      tokens).
- [x] `main_window.py`: `_apply_theme` delega en `theme.apply_theme`
      (conservando QSS de toolbars y refresco de iconos); persistencia de la
      elección de tema en `_apply_preferences`.
- [x] `shell/assembly.py`: carga del tema persistido (reemplaza el
      `"light"` hardcodeado como valor final).

## 4. Persistencia (M21)

- [x] `platform/settings.py`: `UiPreferences` + `load_ui_preferences` /
      `save_ui_preferences` (clave `ui/theme`, normalización a `light`).
- [x] `platform/__init__.py`: export aditivo de la nueva API.

## 5. Infraestructura de iconos (preparación Fase 2)

- [x] `src/chemuson/gui/theme/icon_provider.py`: `IconProvider` SVG→QIcon/
      QPixmap, tinte por `currentColor`, caché por (nombre, color, tamaño),
      HiDPI, glifos dinámicos, fallo visible sin excepción.
- [x] Documentar el contrato de migración de `icons.py` en el provider
      (Fase 2): `icons.py` intacto en esta fase.

## 6. Registro arquitectónico

- [x] `architecture/modules.yml`: M08 paths/`internal_api`/tests de
      `gui/theme`; M21 `public_api` de preferencias UI.

## 7. Tests

- [x] `tests/test_ui_theme_foundation.py`: tokens light/dark, colores
      válidos, métricas, resolución de nombres (incl. `system` y basura),
      `apply_theme` light/dark sin excepciones (QSS + QPalette coherentes),
      compatibilidad de la fachada `styles.py`, ventana real con ambos
      temas, IconProvider (caché/HiDPI/glifos/ausente).
- [x] `tests/test_platform_settings.py`: round-trip y normalización de
      `ui/theme`.
- [x] Smoke Qt offscreen de la ventana real + capturas light/dark en
      `docs/ui-modernization/foundation-shots/`.

## 8. Verificación final

- [x] `openspec validate 2026-09-24-modernize-ui-theme-foundation --strict`.
- [x] Suite completa + targeted tests + compileall + ruff scoped +
      `git diff --check` sin regresiones contra la baseline.
- [x] Revisión del diff: solo archivos del alcance; nada de Clean2D/química.
- [ ] Commit `Add UI theme foundation and design tokens` en
      `ui/modernization`; push si el acceso remoto funciona.
