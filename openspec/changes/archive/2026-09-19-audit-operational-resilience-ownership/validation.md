# Validation

## Automated validation

- Focused operational-resilience audit test: `1 passed`.
- Full architecture suite: `266 passed`.
- Full regression suite: `1492 passed, 55 skipped`.
- Compileall, scoped Ruff and `git diff --check`: passed.
- OpenSpec change validation: passed.
- Full OpenSpec validation: `33 passed, 1 failed` on the pre-existing
  `application-composition-root` requirement without `SHALL`/`MUST`.

## Manual validation

Not applicable: this phase changes no runtime code or GUI behavior.
