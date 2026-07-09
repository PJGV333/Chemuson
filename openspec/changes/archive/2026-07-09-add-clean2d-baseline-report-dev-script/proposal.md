# Proposal: Add Clean 2D Baseline Report Developer Script

## Why
Clean 2D baseline reports and review summaries are available as test-owned helpers, but developers need a simple console entry point to generate, compare, and review reports during local investigation without activating CI gates.

## What Changes
- Add a developer/test-only script for Clean 2D baseline reports.
- Support `write`, `compare`, and `review` commands.
- Serialize reports as JSON-stable output.
- Compare report files with existing report diff rules.
- Review report diffs with the existing baseline diff review policy.
- Return simple exit codes for equivalent, changed, and error outcomes.

## Non-Goals
- No production changes.
- No layout, ranking, selection, UI, backend, RDKit, or CoordGen changes.
- No public document format changes.
- No geometry gates or CI enforcement.
- No geometry repairs.
- No large `engine.py` refactor.

## Success Criteria
- OpenSpec validates before implementation.
- `write` creates a JSON-stable report at the requested path.
- `compare` prints stable diff JSON and returns `0` for equivalent reports, `1` for differences.
- `review` prints stable review summary JSON and returns `0` when no review items exist, `1` when they do.
- Argument, read, write, and JSON errors return `2`.
