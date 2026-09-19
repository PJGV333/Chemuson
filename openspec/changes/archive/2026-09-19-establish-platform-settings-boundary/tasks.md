## Recognition

- [x] Inventory QSettings, settings, preferences, theme and resource uses.
- [x] Classify application configuration, GUI state, document state and runtime state.
- [x] Confirm M21 is the next available catalog ID.

## Tests first

- [x] Add behavior tests for coercion and naming/numbering round-trips.
- [x] Observe the expected RED import failure before implementation.
- [x] Add architecture tests for platform isolation and legacy shims.

## Implementation

- [x] Create M21 platform settings and resource modules.
- [x] Delegate GUI preference policy to M21 while preserving wrappers and keys.
- [x] Preserve the historical utils resource import through a shim.
- [x] Update catalog, docs and dependency contracts.

## Validation

- [x] Run focused tests, architecture tests and full regressions.
- [x] Run compileall, Ruff, OpenSpec and diff checks.
