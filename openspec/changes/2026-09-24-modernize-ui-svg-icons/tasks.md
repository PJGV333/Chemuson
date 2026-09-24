# Tasks: Sistema de iconos SVG de producción

## 0. Preparación (OpenSpec + baseline)

- [x] Crear el cambio OpenSpec `2026-09-24-modernize-ui-svg-icons`
  (proposal, design con inventario 1:1, tasks, spec).
- [x] `openspec validate 2026-09-24-modernize-ui-svg-icons --strict`
  (antes de implementar).
- [x] Baseline en `baseline.md` (git status, compileall, pytest
  collect-only, pytest completo, ruff scoped) con el entorno reproducible
  de la Fase 1.

## 1. Set de SVG estáticos (55)

- [x] Crear `src/chemuson/gui/theme/icons/` con los 24 iconos genéricos
  (pointer, eraser, pan, rotate-left, rotate-right, flip-horizontal,
  flip-vertical, zoom-in, zoom-out, chain, lasso, corner, frame,
  rounded-frame, tlc, electrophoresis, doc-new, doc-open, doc-save, undo,
  redo, copy, paste, clean) en lenguaje 24×24/1.75/round/currentColor.
- [x] Crear los 11 `bond-*` y los 17 `arrow-*`.
- [x] Crear `energy-levels` (pastels bakes), `molecular-orbital`,
  `wavy-anchor`.
- [x] `LICENSE.txt` con origen/atribución del set.

## 2. Generadores dinámicos

- [x] `src/chemuson/gui/theme/icon_svg.py`: builders puros `atom`,
  `sphere`, `charge`, `electrons`, `radical`, `ring`, `ring-template`,
  `energy-boxes` (SVG 24×24, `currentColor` salvo colores de dominio) +
  registro `BUILDERS`.
- [x] `icon_provider.py`: `icon_dynamic` / `pixmap_dynamic` (mismo
  render/tinte/caché/HiDPI; clave `dyn:<key>:<params>`), sin romper la API
  de Fase 1.

## 3. Fachada `gui/icons.py`

- [x] Reescribir `icons.py` como fachada 1:1 sobre el provider (mismos
  nombres, firmas, `ICON_SIZE`, `ATOM_COLORS`, `set_icon_theme`,
  `icon_*_color` → tokens).
- [x] Eliminar el pixel-loop de tint y `QIcon.fromTheme` (undo/redo → SVG).
- [x] Verificar que los 4 callers (`toolbar.py`, `text_toolbar.py`,
  `main_window.py`, `main_window_ui_builder.py`) no necesitan cambios.

## 4. Registro y tests

- [x] `architecture/modules.yml`: M08 registra `theme/icons/` y
  `theme/icon_svg.py`.
- [x] `tests/test_ui_svg_icons.py`: inventario/paridad, validez de SVG,
  provider (tamaños, HiDPI, caché, temas, fallo seguro), glifos
  dinámicos, compatibilidad de la fachada, smoke light/dark de la ventana.
- [x] Tests targeted verdes: `pytest tests/test_ui_svg_icons.py
  tests/test_ui_theme_foundation.py -q`.

## 5. Verificación final

- [x] Suite completa: mismos resultados que baseline (sin regresiones
  nuevas; los 4 fallos preexistentes permanecen sin tocar).
- [x] `python -m compileall src tests tools packaging` sin errores.
- [x] `ruff check src tests tools packaging --select
  F401,F811,F821,E722,E741` sin findings nuevos.
- [x] `git diff --check` limpio.
- [x] Smoke Qt offscreen: ventana real con ambos temas; capturas light/dark
  en `docs/ui-modernization/icon-phase-shots/` + README de comparación.
- [x] Commit limpio exclusivo de la Fase 2; push si la autenticación
  permite (si no, dejar el commit local listo y reportarlo).
