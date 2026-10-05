## 1. OpenSpec and baseline

- [x] 1.1 Capture the clean Phase 2 commit, compileall, test collection, scoped fast tests, Ruff, and known full-suite time limit before changes.
- [x] 1.2 Validate this proposal, design, and evaluator requirements with strict OpenSpec validation.

## 2. Non-mutating evaluator

- [x] 2.1 Add a `tools/` report builder that gates on successful isolated-ChemIO M23 results and never calls Clean2D for failures.
- [x] 2.2 Measure initial and selected-candidate geometry with existing M02 APIs and preserve engine state/reason without using metrics as gates.
- [x] 2.3 Add an explicit CLI requiring endpoint/model/description; emit finite, JSON-safe reports without prompt, credential, or raw provider error data.
- [x] 2.4 Confirm evaluation leaves the validated graph and its coordinates unchanged and has no canvas/document integration.

## 3. Verification and boundaries

- [x] 3.1 Add offline fake-provider/engine tests for success/failure routing, report shape, privacy, and serialization.
- [x] 3.2 Add one short, bounded test through the real Clean2D engine on a small deterministic graph.
- [x] 3.3 Verify M02 remains free of M23 imports and no new runtime module/dependency is introduced.
- [x] 3.4 Run focused tests under a five-minute shell timeout, compileall, scoped Ruff, strict OpenSpec validation, and diff checks; record results.
- [x] 3.5 Do not rerun full pytest: the known full-suite baseline takes 19:26 and exceeds the requested maximum test duration.
