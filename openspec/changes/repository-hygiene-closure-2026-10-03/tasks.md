# Tasks: Repository hygiene and historical consolidation

## 1. Baseline and campaign memory
- [x] Record baseline measurements, test results, HEAD and status.
- [x] Audit Git history, OpenSpec archives, prior reports and all remote branches.
- [x] Write `docs/history/CAMPAIGNS.md` with campaign outcomes, decisions, lessons and branch SHAs.
- [x] Write `docs/history/REPOSITORY_POLICY.md`.

## 2. Conservative audits
- [x] Audit all remote/local branches; classify each and record safe deletion candidates and unique work.
- [x] Audit `src/sys` references and package/test/asset loading; remove only if unused.
- [x] Audit code/API candidates with imports/AST and dynamic consumer checks; remove only SAFE_DELETE candidates.
- [x] Audit UI/OpenSpec screenshots, `tests/archive`, historical assets and redundant reports; preserve Markdown/canonical specs and useful evidence.

## 3. Cleanup and validation
- [x] Remove confirmed artifacts/files and add only precise ignore rules.
- [x] Run compileall, architecture tests, targeted UI tests, full suite, scoped Ruff, OpenSpec strict, diff check and Qt offscreen smoke; investigate any new failure and stop if regression/asset loss.
- [x] Re-measure tree, source/docs/tests/OpenSpec, `.git`, pack size, bytes/files and blob inventory; record future history-rewrite estimate without executing it.

## 4. Branch lifecycle and closure
- [x] Push `maintenance/repository-hygiene-closure` normally.
- [x] Delete only documented remote/local SAFE_TO_DELETE branches that are ancestors of main; preserve protected and all unique/active/unknown branches.
- [x] Fetch with prune and record resulting remote inventory.
- [x] Write `docs/history/REPOSITORY_CLEANUP_2026-10-03.md` and complete this checklist only after all gates pass.
- [x] Do not merge to main, rewrite history, or begin another product phase.
