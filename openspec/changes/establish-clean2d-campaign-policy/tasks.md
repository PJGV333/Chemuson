# Tasks: Establish Clean2D Campaign Policy

## 1. Repository and baseline

- [x] Fetch `origin`, verify a clean starting tree, record branch and base commit.
- [x] Create the dedicated branch `clean2d/campaign-policy`.
- [x] Read the active and archived Clean2D OpenSpecs, architecture documents, production Clean2D modules, and regression corpus.
- [x] Capture the required baseline commands in `baseline.md`.

## 2. Master OpenSpec

- [x] Create `proposal.md`, `design.md`, `tasks.md`, and the campaign-policy spec.
- [x] Reuse existing result-state, rejection-reason, diagnostic, snapshot, corpus, metric, and baseline vocabulary.
- [x] Define the simple/medium/large/complex-scale policy and preserve-only/no-op rule.
- [x] Define hard constraints, hard gates, soft metric vector, taxonomy, case identity, baseline workflow, determinism, observability, performance, and external-backend policy.
- [x] Define all nine sequential campaigns and require one OpenSpec per future campaign.
- [x] Define promotion gates, rollback, experimental routing, visual review, and production-ready criteria.

## 3. Roadmap and contract guard

- [x] Add `docs/clean2d/CAMPAIGN.md` with the required human roadmap sections.
- [x] Add `tests/architecture/test_clean2d_campaign_policy.py` for structural and normative policy checks; verify production scope by diff review rather than Git-history state.
- [x] Keep the change documentary/contractual; do not modify production Clean2D or existing specifications.

## 4. Validation

- [x] Run `openspec validate establish-clean2d-campaign-policy --strict` (valid).
- [x] Run `openspec validate --all --strict` (executed; repository has pre-existing placeholder-purpose failures in other specs).
- [x] Run the new contract test and the existing Clean2D suite (1,536 passed, 20 skipped, 4 failures in untouched existing files: candidate generation and stereo import).
- [x] Run compileall, targeted Ruff, `git diff --check`, and final status.

## 5. Commit and publication

- [x] Create one commit named `Establish Clean2D campaign policy`.
- [x] Publish the branch if the configured Git remote permits it.
- [ ] Do not archive or merge this OpenSpec.
