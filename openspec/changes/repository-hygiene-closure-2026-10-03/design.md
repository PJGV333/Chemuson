# Design: Repository hygiene closure

## Baseline

The immutable start is `1db4f63b52af79247745b3a8a220fb728348218c` on clean `main`, equal to `origin/main`; work is on `maintenance/repository-hygiene-closure`. Exact measurements and the most recent full-suite result for this exact source tree are in `baseline.md`; detailed object and branch inventories are under `/tmp/chemuson-*`.

## Procedure

1. Preserve current canonical requirements and read the repository instructions. Record disk, tracked-tree, Git pack/object, file, PNG, blob and branch baselines without decoding binary blobs as UTF-8.
2. Summarize campaigns from committed OpenSpec, reports and branch diffs before any historical branch/artifact deletion. Distinguish merged, active, superseded, unique and unknown work using ancestry and patch/content evidence, not names.
3. Build a conservative candidate table: SAFE_DELETE / KEEP / NEEDS_REVIEW with evidence for each decision. Specifically verify `src/sys` through imports, textual references, package manifests, assets, tests and entry points. Use a path-specific ignore rule only if removal is justified.
4. Keep canonical `openspec/specs/`; retain Markdown archives. Prune archived binary evidence only when it is redundant with final canonical evidence or its technical outcome is adequately described elsewhere. Never delete active behavior tests merely to reduce size.
5. Implement only confirmed removals and documentation. Avoid product-code refactors; if any code candidate is ambiguous, preserve and report it.
6. Re-measure and run compileall, architecture, targeted UI, full suite, Ruff, all strict OpenSpec, `git diff --check`, and Qt offscreen smoke. Compare failures with the baseline; stop on new functional failures or asset loss.
7. Push the maintenance branch normally. Only then delete remote branches documented as fully merged ancestors. Preserve active/unique/unknown branches and list recommendations for owner review. Fetch with prune and measure the final branch inventory.
8. Do not rewrite Git history. Estimate the potential savings and affected refs for a future, separately approved filter-repo operation.

## Accounting

Report both tracked logical file bytes and directory disk usage, and distinguish `.git` total size from packed objects. The ignored local `.venv` is environment state, not a tracked deliverable; measure and disclose it separately rather than deleting it opportunistically. Keep detailed lists in `/tmp`, with only summaries/top 30 in the report.
