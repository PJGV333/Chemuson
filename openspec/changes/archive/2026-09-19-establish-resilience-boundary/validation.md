# Validation

## Automated validation

- Focused resilience/startup tests: `55 passed`.
- Full architecture suite: `263 passed`.
- Full regression suite: `1489 passed, 55 skipped`.
- `python -m compileall -q src tests tools packaging`: passed.
- Scoped Ruff for M22, shims and consumers: `All checks passed!`.
- OpenSpec change validation: passed.
- `git diff --check`: passed.

## Manual validation

- `QT_QPA_PLATFORM=offscreen uv run --no-sync python` launched
  `ChemusonWindow`, ran the Qt event loop, and exited with
  `qt_resilience_smoke_exit=0`.
- The autosave subprocess isolation contract confirmed no eager PyQt6, GUI,
  ChemIO or RDKit imports through the historical shim.
