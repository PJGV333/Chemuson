# Validation

## Automated validation

- Focused clipboard policy and codec regressions passed.
- Related copy/paste and image clipboard regressions passed.
- Full regression baseline passed: 1477 passed, 55 skipped.
- Compileall passed.
- Targeted Ruff passed; the repository-wide required selection retains the pre-existing Clean2D F401 documented in `AGENT_REPORT.md`.
- The archived OpenSpec was validated successfully before archival.

## Manual validation

Completed: Qt smoke testing exercised copy/paste and clipboard behavior. No
behavioral regression was observed.
