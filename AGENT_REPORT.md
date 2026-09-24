# AGENT_REPORT

## Scope

Phase 15 closed the selection-helper extraction and narrowed M20 to the five
canonical helpers under `gui/editor2d/selection/`. The parent
`gui/editor2d/__init__.py` remains a logic-free namespace owned by M08, and
M09 retains import-only compatibility shims under `gui/canvas`.

## Closure

The following OpenSpecs were archived on 2026-09-19 after the reported manual
Qt smoke test completed:

- `2026-09-19-extract-canvas-selection-hit-testing`
- `2026-09-19-extract-canvas-selection-overlays-and-handles`
- `2026-09-19-extract-canvas-selection-clipboard-policy`
- `2026-09-19-move-canvas-selection-helpers-to-editor2d`

Manual validation exercised selection, overlay handles, copy/paste and legacy
import consumers. The archived task and validation records mark all work
complete.

## Baseline and validation deviations

- Baseline pytest collection: 1532 tests collected.
- Latest complete pytest run passes: 1477 passed, 55 skipped.
- Compileall passes for `src tests tools packaging`.
- Targeted checks for the migrated helpers, shims, consumers and architecture
  tests pass.
- The repository-wide required Ruff selection reports exactly one pre-existing
  out-of-scope F401 in
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`).
  It was not modified because this phase does not cover Clean2D tests.
- `openspec validate --all --strict` retains exactly one pre-existing unrelated
  failure in `spec/application-composition-root`: its first requirement lacks a
  SHALL or MUST keyword. No exception or catalog workaround was added.

Neither deviation is caused by the selection ownership or extraction work.

## Remaining module-boundary wave

The descendant branch `architecture/remaining-module-boundaries` completed the
following audited boundaries:

- Phase A: M05/M06/M07/M14/M16/M17/M18/M19/M20 audited; cohesive modules kept intact.
- Phase B: M21 `platform.settings` owns application preferences and packaged resources.
- M22 `resilience` owns canonical autosave, recovery filesystem policy and crash logging; historical utils paths remain import-only shims.
- Phase D: M14 update subsystem audited; M23 remains reserved.
- Phase E: M19 composition root audited; M24 remains reserved.
- Phase F: operational resilience ownership documented across M22, M08/M10, M14 and M19.

Final branch validation:

- Architecture suite: `266 passed`.
- Full regression suite: `1492 passed, 55 skipped`.
- Compileall: passed.
- Scoped Ruff for changed implementation/tests: passed.
- Final Qt offscreen smoke: `qt_resilience_smoke_exit=0`.
- Full OpenSpec validation: `33 passed, 1 failed`, the same pre-existing
  `application-composition-root` requirement issue.
- Full required Ruff selection: the same pre-existing Clean2D `math` F401.

All phase OpenSpecs are archived and the worktree is clean after the final
operational-resilience audit commit.

## Final boundary reconciliation

- M15 is documented and tested as compatibility shims only.
- M19 current and target dependencies are explicitly M08/M18/M22.
- M20 owns only `gui/editor2d/selection/`; M09 owns legacy canvas shims.
- M21 lists `platform/__init__.py`, `settings.py` and `resources.py` explicitly.
- M22 now owns GUI-free recovery filesystem policy in `resilience/recovery.py`.
- Final focused recovery/platform tests: `20 passed`.
- Final architecture suite: `268 passed`.
- Final full regression suite: `1496 passed, 55 skipped`.
- Compileall and scoped Ruff passed; global OpenSpec validation is `35 passed, 0 failed`.
- Qt offscreen smoke exited with `qt_recovery_smoke_exit=0`.

## UI modernization — Fase 1: theme foundation (2026-09-24, branch ui/modernization)

OpenSpec: `2026-09-24-modernize-ui-theme-foundation` (validado `--strict`).
Implementación de la Fase 0+1 de `docs/ui-modernization/PLAN.md` sobre la
referencia visual aprobada del spike PyQt6 (commit `59e977d`).

### Decisiones y desviaciones registradas

- **Tokens: el spike es la referencia, no la tabla textual de PLAN.md §2.2**
  (diferencias de valores y de vocabulario; registrado en el proposal §6 y
  design D1). Único ajuste sobre el spike: los rellenos suaves del tema
  oscuro pasaron de `rgba(r,g,b,a)` a hex `#AARRGGBB` (color idéntico)
  porque `QColor` no parsea `rgba()` (mantiene "todo token es color válido";
  QSS acepta ambos formatos).
- **`styles.py` como fachada**: PLAN.md Fase 1 hablaba de "función única";
  en producción existen 2 generadores + 2 constantes con callers reales
  (`main_window.py`, `toolbar.py`); la fachada conserva los 4 nombres y las
  paletas legadas como alias de tokens (deprecation documentada).
- **Persistencia del tema**: antes no existía (`assembly.py` hardcodeaba
  `"light"`). Se añade `ui/theme` vía M21 (`UiPreferences` + load/save);
  `"system"` queda aceptado/normalizado por la API y persistible, pero sin
  opción visible en Preferencias (viene con la fase de ajustes).
- **IconProvider**: implementado como infraestructura (Fase 2 migra
  `icons.py`); `icons.py` intacto según instrucción.
- **Sin nuevas dependencias**: solo PyQt6 (QtSvg ya disponible).

### Entorno de tests (quirk)

- Worktree sin `.venv`; el venv del checkout principal tiene el stack
  runtime pero no pytest/ruff. Validación con entorno efímero:
  `uv run --no-project --offline --python <venv principal>/bin/python --with
  pytest --with ruff --with PyQt6 --with numpy --with Pillow --with rdkit
  --with certifi --with PyYAML` (todo resuelto desde la caché de uv; sin
  instalaciones nuevas). Con `HOME` aislada es necesario fijar también
  `UV_CACHE_DIR` (uv usa `$HOME/.cache/uv`).
- Baseline (ver `openspec/changes/2026-09-24-modernize-ui-theme-foundation/baseline.md`):
  1583 passed / 20 skipped / 4 failed preexistentes (candidate generation +
  stereo import) y 1 ruff F401 preexistente en un test de Clean2D (fuera de
  alcance; no corregido por política de no refactors oportunistas).

### Verificación

- `openspec validate 2026-09-24-modernize-ui-theme-foundation --strict`: OK.
- Tests de la fundación: `tests/test_ui_theme_foundation.py` (tokens,
  resolución light/dark/system, aplicación sin excepciones, QSS/QPalette,
  fachada `styles.py`, ventana real, IconProvider caché/HiDPI) +
  `tests/test_platform_settings.py` (round-trip `ui/theme`).
- Smoke Qt offscreen de la ventana real con ambos temas:
  `SMOKE: OK` (capturas en `docs/ui-modernization/foundation-shots/`).
- Ámbito del diff restringido a: `gui/theme/` (nuevo), `styles.py`,
  `main_window.py`, `shell/assembly.py`, `platform/`, `architecture/modules.yml`,
  tests y OpenSpec/docs. Intactos: `clean2d/`, `chemname/`,
  `chemio/persistence.py`, `gui/canvas/`, `gui/icons.py`, `gui/toolbar.py`,
  `gui/docks.py`, `gui/style.py`, `gui/dialogs/`.
