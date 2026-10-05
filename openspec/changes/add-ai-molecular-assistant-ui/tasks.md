## 1. OpenSpec and baseline

- [x] 1.1 Record the clean branch/HEAD and exact Phase 2 baseline outputs in this change.
- [x] 1.2 Validate the proposal, design, and UI requirements with strict OpenSpec validation.

## 2. GUI entry point and request flow

- [x] 2.1 Add one AI structure-generation QAction to the Structure menu and register it in the existing command palette.
- [x] 2.2 Add a minimal modeless request/review UI with a description, explicit endpoint/model, optional masked API key, and transient configuration only.
- [x] 2.3 Add a controller/worker path that runs M23 generation and isolated ChemIO validation off the GUI thread; suppress abandoned/late results and retain bounded timeout behavior.
- [x] 2.4 Present controlled failures by stable status/reason only; never expose raw provider diagnostics or insert a failed result.
- [x] 2.5 Show a successful result's provenance, exact SMILES, and semantic caveat; require explicit approval before insertion through the existing canvas undo macro.

## 3. Focused verification and architecture

- [x] 3.1 Add offline UI/controller tests for discovery, validation/config rejection, async delivery, result preview, and transient secret handling.
- [x] 3.2 Add editor-state snapshots proving failure, decline, and abandoned/late results preserve graph, selection, undo index, and dirty state.
- [x] 3.3 Verify approved insertion is selected, represented as one undo step, and restored by Undo/Redo.
- [x] 3.4 Update `architecture/modules.yml` and architecture tests to record M10's one-way M23 dependency while preserving M23's M00/M01-only dependencies and M02's AI independence.
- [x] 3.5 Run focused UI/controller/canvas tests, architecture tests, compileall, collection, scoped Ruff, strict OpenSpec validation, and `git diff --check`; record baseline-relative results and known findings.
- [ ] 3.6 Full `pytest -q` was not rerun: the captured baseline takes 19:26 (over the requested 5–10 minute test cap) and has one unrelated pre-existing CompChem failure; this gate is explicitly deferred for this phase.
