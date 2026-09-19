# Validation

## Automated validation

- Focused hit-testing regressions passed.
- Combined selection, semantic-diagram and tightness regressions passed: 29 tests.
- Full regression baseline passed: 1477 passed, 55 skipped.
- Compileall passed.
- Targeted Ruff passed; the repository-wide required selection retains the pre-existing Clean2D F401 documented in `AGENT_REPORT.md`.
- The archived OpenSpec was validated successfully before archival.

## Manual validation

Completed: Qt smoke testing exercised selection hit testing and legacy import
consumers. No behavioral regression was observed.
