# Validation

## Automated validation

- Architecture audit contract: `2 passed`.
- Full architecture suite: `256 passed`.
- Full regression suite: `1479 passed, 55 skipped`.
- `python -m compileall -q src tests tools packaging`: passed.
- Scoped Ruff for the added architecture test: passed.
- `openspec validate audit-existing-domain-module-boundaries --strict`: valid.
- Repository-wide OpenSpec validation retains the known unrelated
  `application-composition-root` requirement failure documented in
  `AGENT_REPORT.md`.
- `git diff --check`: passed.

## Manual validation

Not applicable: this phase changes no Qt behavior or production runtime code.
