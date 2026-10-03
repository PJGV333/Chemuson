# Proposal: Repository hygiene and historical consolidation

## Purpose

Leave ChemUSON maintainable after its architecture, Clean2D and UI campaigns: preserve actionable technical history in a concise document, remove confirmed accidental/generated artifacts, and retire only remote branches proven safe to delete.

## Scope

- Capture before/after repository, Git-object, branch and test metrics.
- Create `docs/history/CAMPAIGNS.md`, `REPOSITORY_POLICY.md`, and a dated cleanup report.
- Audit every remote branch by ancestry and content; document every branch considered for deletion before deleting it.
- Audit dead code, `src/sys`, archive screenshots, tests/archive, and historical documents conservatively. Delete only candidates proven unused/redundant; retain canonical specs and current tests.
- Run compile, architecture, UI, full-suite, Ruff, OpenSpec, diff and Qt smoke checks.
- Push only `maintenance/repository-hygiene-closure`; do not merge it to `main`.

## Non-goals

No product feature, chemical behavior, Clean2D algorithm, canvas/editor, CMSN, template chemistry, export, public API, version, Git history rewrite, force-push or automatic merge. Do not delete `main`, `gh-pages`, active Clean2D branches, or any branch with unique work that has not been reviewed and documented.

## Deletion gates

A file is removable only after consumer checks (imports, AST, dynamic hooks, signals/actions, packaging/assets, tests and documentation). A branch is removable only after its HEAD is proven an ancestor of `origin/main`, the branch SHA and absorbed campaign are documented, the maintenance branch is pushed, and validation is green. Unique non-ancestor work is preserved and receives an explicit recommendation.
