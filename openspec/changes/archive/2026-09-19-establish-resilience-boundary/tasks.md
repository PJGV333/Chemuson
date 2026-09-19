## Recognition

- [x] Audit crash logging, autosave, recovery and failure-isolation consumers.
- [x] Confirm autosave and crash reporting form a cohesive generic boundary.

## Tests first

- [x] Add canonical-path and import-isolation architecture tests.
- [x] Observe the expected RED failure before moving canonical files.
- [x] Preserve existing autosave and startup tests through shims.

## Implementation

- [x] Move canonical autosave and crash reporting into M22.
- [x] Keep import-only compatibility shims under M15.
- [x] Update bootstrap, tab manager and recovery imports.
- [x] Update catalog, docs and module-count contracts.

## Validation

- [x] Run focused, architecture and full regression suites.
- [x] Run compileall, scoped Ruff, OpenSpec, diff and Qt smoke checks.
