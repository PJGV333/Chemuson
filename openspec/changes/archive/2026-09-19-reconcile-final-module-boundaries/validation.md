# Validation

## Automated validation

- Recovery/platform/composition focus: `20 passed`.
- Full architecture suite: `268 passed`.
- Full regression suite: `1496 passed, 55 skipped`.
- `python -m compileall src tests tools packaging`: passed.
- Scoped Ruff for changed implementation/tests: `All checks passed!`.
- `openspec validate reconcile-final-module-boundaries --strict`: valid.
- `openspec validate --all --strict`: `35 passed, 0 failed`.
- `git diff --check`: passed.

## Manual validation

- `QT_QPA_PLATFORM=offscreen uv run --no-sync python` launched
  `ChemusonWindow`, ran the Qt event loop, and exited with
  `qt_recovery_smoke_exit=0`.
- The Qt platform plugin emitted the non-blocking warning
  `This plugin does not support propagateSizeHints()`.
- The repository-wide required Ruff selection still reports exactly one
  pre-existing out-of-scope F401 in
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`).
- No credentials or secrets were stored; sensitive values are represented as
  `[REDACTED]`.
