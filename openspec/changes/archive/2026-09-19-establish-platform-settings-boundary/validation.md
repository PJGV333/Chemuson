# Validation

## Automated validation

- Focused platform settings/resource tests: `52 passed`.
- Full architecture suite: `259 passed`.
- Full regression suite: `1485 passed, 55 skipped`.
- `python -m compileall -q src tests tools packaging`: passed.
- Scoped Ruff for M21 consumers and tests: `All checks passed!`.
- OpenSpec change validation: passed.
- Full OpenSpec validation: `29 passed, 1 failed` on the pre-existing
  `application-composition-root` requirement without `SHALL`/`MUST`.
- `git diff --check`: passed.

## Manual validation

- `QT_QPA_PLATFORM=offscreen uv run python` launched `ChemusonWindow`, ran the
  Qt event loop, and exited with `qt_smoke_exit=0`.
- Preference round-trips, resource imports and legacy resource shim tests passed.
