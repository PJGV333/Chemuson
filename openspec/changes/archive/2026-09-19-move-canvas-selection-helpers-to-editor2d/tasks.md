## Recognition

- [x] Inspect all helper definitions, imports, consumers and catalog ownership.

## Tests first

- [x] Add canonical package, shim and catalog ownership tests.
- [x] Verify the tests fail before moving implementations.

## Migration

- [x] Create M20 `gui/editor2d/selection/` and move the five implementations.
- [x] Replace old modules with compatibility shims.
- [x] Update M09 imports to canonical M20 paths.

## Architecture and docs

- [x] Register M20 and remove helper ownership from M09.
- [x] Update M09 and add M20 documentation.
- [x] Record manual GUI validation status.

## Validation

- [x] Run focal and full regression suites.
- [x] Run compileall, Ruff, OpenSpec and diff checks.
- [x] Commit and archive the migration after GUI validation.
