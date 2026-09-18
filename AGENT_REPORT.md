# AGENT_REPORT

## Scope

Phase 15 moved the five extracted canvas selection helpers to M20
(`gui/editor2d`) and retained import-only shims under M09 (`gui/canvas`).

## Validation deviations

- The complete pytest suite passes: 1474 passed, 55 skipped.
- Changed-file compileall/Ruff checks pass. The repository-wide required Ruff
  selection reports one pre-existing out-of-scope F401 in
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`).
  It was not modified because this phase does not cover Clean2D tests.
- `openspec validate --all --strict` reports one pre-existing unrelated failure
  in `spec/application-composition-root`: its first requirement lacks a SHALL
  or MUST keyword. The phase-specific OpenSpec validates successfully, as do
  the other 27 cataloged items.

Neither deviation is caused by the M20 migration. No architecture exception or
catalog workaround was added. Manual Qt GUI validation remains pending before
archiving the migration change.
