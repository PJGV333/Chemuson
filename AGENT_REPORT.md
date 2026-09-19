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
- Phase C: M22 `resilience` owns autosave and crash logging; historical utils shims remain.
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
