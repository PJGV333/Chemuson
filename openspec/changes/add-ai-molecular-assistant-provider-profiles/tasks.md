## 1. OpenSpec and baseline

- [x] 1.1 Capture the clean Phase 3 commit, compileall, collection, focused UI/provider/architecture tests, Ruff, and full-suite time limit before changes.
- [x] 1.2 Validate the proposal, design, and provider profile/UI deltas with strict OpenSpec validation.

## 2. Provider profile catalog and UI

- [x] 2.1 Add immutable M23 metadata for custom, OpenAI, LM Studio, and llama.cpp-compatible endpoint profiles; preserve generic config defaults.
- [x] 2.2 Add a profile selector to the existing dialog with editable endpoint/model fields, provider-aware key labels, and key clearing when profiles change.
- [x] 2.3 Pass the selected profile ID through the existing controller config into M23 provider provenance.
- [x] 2.4 Enforce the OpenAI key requirement before starting a worker; keep local/custom credentials optional and transient.

## 3. Verification and boundaries

- [x] 3.1 Add offline profile and fake-transport tests for URLs, model IDs, auth headers, and existing strict decoder behavior.
- [x] 3.2 Add dialog/controller tests for profile defaults, editable model IDs, key requirement, key clearing, and preserved asynchronous review/insertion.
- [x] 3.3 Verify the M23 catalog/API entry and M10→M23 direction; add no runtime dependency or native provider protocol.
- [x] 3.4 Run focused tests under a five-minute timeout, compileall, test collection, scoped Ruff, strict OpenSpec validation, and diff checks; document global baseline findings.
- [x] 3.5 Do not run provider/network or manual model-server tests. Do not rerun full pytest because its recorded duration is 19:26.
