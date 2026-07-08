# Tasks: Create Clean 2D Regression Corpus

## 1. OpenSpec Artifacts
- [x] Update `proposal.md` to state that this change creates regression infrastructure, not geometry improvements.
- [x] Update `design.md` with explicit layers A-F.
- [x] Update `tasks.md` with small implementation phases.
- [x] Add `specs/clean-2d-regression-corpus/spec.md` with verifiable requirements.
- [x] Run `openspec validate create-clean2d-regression-corpus --strict` before implementation.

## 2. Corpus Registry
- [x] Add a test-only regression case dataclass.
- [x] Add a stable registry with unique case names.
- [x] Require each case to declare family, mode, target, expected contract states, and known-delicate metadata.

## 3. Initial Molecule Cases
- [x] Add a simple acyclic valid graph case.
- [x] Add a simple aromatic valid graph case.
- [x] Add a fused-ring or fused-aromatic representative case without improving its geometry.
- [x] Add a partial-selection or boundary-sensitive case.
- [x] Add a dirty-coordinate valid graph case.

## 4. Contract Assertion Helpers
- [x] Add helpers that capture chemical identity before and after Clean 2D.
- [x] Assert atom ids, bond ids, elements, endpoints, bond orders, and stable stereo metadata are preserved.
- [x] Assert Clean 2D execution does not raise uncontrolled exceptions.
- [x] Assert result state is one of the case's declared expected states.

## 5. Snapshot Serialization Tests
- [x] Add debug snapshot generation for case metadata, result contract, identity signatures, and optional diagnostics.
- [x] Assert every snapshot is JSON-serializable.
- [x] Keep snapshots test-owned and separate from public document format.

## 6. Regression Corpus Tests
- [x] Add pytest coverage that iterates all registered cases.
- [x] Add explicit tests for unique case names and required case metadata.
- [x] Add coverage for known-delicate handling without forcing geometry fixes.

## 7. Validation
- [x] Run `openspec validate create-clean2d-regression-corpus --strict` after implementation.
- [x] Run the new regression corpus pytest tests.
- [x] Run relevant existing Clean 2D tests when feasible.

## 8. Follow-Up OpenSpec Roadmap
- [x] Document recommended follow-up changes for expanding the corpus, defining metric contracts, and improving geometry based on measured baselines.
