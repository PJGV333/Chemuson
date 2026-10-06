# Design

## Decisions

1. `pyproject.toml` remains the canonical runtime dependency declaration. `requirements.txt` mirrors it, including Pillow; `requirements-dev.txt` contains only pytest, Ruff, and PyYAML. README documents editable installation and an arbitrary venv name.
2. The evaluator uses the existing `summarize_clean2d_candidates()` contract. Summary values are normalized through a strict JSON-safe conversion; no ranking, selection, or Clean2D implementation is changed.
3. `--api-key-env` has no implicit environment-variable default. The CLI reads a secret only when a variable name is explicitly supplied; the report and diagnostics never contain the key.
4. Existing baseline summaries/validation records are retained. Duplicated full command logs are not versioned; command output stays under `/tmp`.
5. Test coverage uses collection plus bounded targeted suites/shards. The recorded 19:26 baseline makes a monolithic full-suite run incompatible with the current hard timeout policy, so no such run is attempted.

## Verification

Run compileall, collect-only, focused integration tests, architecture tests, strict OpenSpec, scoped Ruff, and `git diff --check`, all under explicit timeouts. Preserve pre-existing Ruff and test findings without suppressing or relabeling them.
