## 1. OpenSpec and baseline

- [x] 1.1 Read active assistant OpenSpecs and capture branch/HEAD, status, compileall, collection, scoped Ruff, and bounded relevant-test baseline.
- [x] 1.2 Record why the monolithic suite is not run: prior runtime 19:26 and current operator cap of 10 minutes.

## 2. Integration hygiene

- [x] 2.1 Remove redundant tracked `baseline-output.log` files from UI, Clean2D-evaluator, and provider-profile changes after confirming baseline.md/validation.md carry the useful results.
- [x] 2.2 Remove the obsolete Phase 5 decision-block from `AGENT_REPORT.md` and preserve concise history in `docs/history/CAMPAIGNS.md`.
- [x] 2.3 Align `requirements.txt`, add `requirements-dev.txt`, and document new-checkout editable installation in README.
- [x] 2.4 Add safe selected/candidate/rejected Clean2D summaries to the evaluator without touching `src/chemuson/clean2d/`.
- [x] 2.5 Make `--api-key-env` opt-in and test that a present `OPENAI_API_KEY` is not used absent the flag.

## 3. Verification and closeout

- [x] 3.1 Run bounded focused tests, architecture, compileall, collect-only, Ruff, strict OpenSpec, and diff check; log each result.
- [x] 3.2 Commit and push only to `origin/ai/molecular-assistant-foundation`; do not merge/rebase/force-push.
