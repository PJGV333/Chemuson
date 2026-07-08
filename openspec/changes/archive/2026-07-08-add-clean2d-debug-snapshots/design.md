## Context

Clean 2D has recently gained a stable decision contract and unified quality reporting. That makes final status and diagnostic data more consistent, but it does not yet provide a single artifact that captures enough information to reproduce difficult Clean 2D failures.

The target use case is debugging complex algorithmic issues in tests and developer runs. Snapshot capture must be purely observational: it can inspect inputs, candidate information, diagnostics, and final state, but it must not influence candidate generation, ranking, selection, geometry normalization, RDKit/CoordGen behavior, MolGraph state, or canvas behavior.

Likely affected areas are the Clean 2D orchestration layer, the quality diagnostic objects already produced by candidate evaluation, and test utilities. Implementation should avoid a global `engine.py` refactor and should prefer a small serialization boundary around data that Clean 2D already computes.

## Goals / Non-Goals

**Goals:**

- Provide an opt-in debug snapshot system for Clean 2D runs.
- Define a stable JSON object that can be serialized, read back, and compared in tests.
- Capture the molecule topology, selection context, initial/final geometry, candidate provenance, candidate diagnostics, and final decision state.
- Support enabling snapshots through an environment variable, an explicit debug parameter, or test-only helper.
- Provide test helpers for writing and reading snapshots.
- Add tests that verify snapshot content and no-op behavior when disabled.

**Non-Goals:**

- Do not implement or change any Clean 2D algorithm.
- Do not change candidate ranking, candidate selection, or fallback behavior.
- Do not change expected geometry.
- Do not refactor `engine.py` globally.
- Do not change RDKit, CoordGen, MolGraph, or canvas code except for minimal integration points needed to observe existing data.
- Do not fix the known algorithmic failures in this change.
- Do not enable snapshots by default.

## Decisions

### Snapshot Format

Use a stable JSON object with explicit `schema` and `version` fields. The object should contain only JSON-serializable values: strings, numbers, booleans, arrays, objects, and null. Coordinates should be represented as keyed atom-coordinate records rather than Python or Qt/RDKit objects.

Rationale: JSON is easy to inspect, store as a fixture, diff across commits, and read from tests without project-specific binary tooling.

Alternative considered: pickle or Python object dumps. Rejected because they are not stable across code changes and are harder to review.

### Observational Capture

Snapshot capture should read data already present during a run. It must not request extra candidates, reorder candidates, recompute geometry for debugging purposes, or call fallback systems that normal execution would not call.

Rationale: The snapshot is for reproducing bugs, so enabling it must not hide or introduce the bug being investigated.

Alternative considered: recomputing diagnostics after the run. Rejected because it risks changing execution timing, accessing different state, or producing diagnostics that were not part of the original decision.

### Opt-In Activation

Support three activation paths:

- Environment variable for developer reproduction runs.
- Explicit debug parameter for direct API or engine-level tests.
- Test-only helper that scopes snapshot capture to a fixture directory.

Rationale: Different debugging workflows need different ergonomics, while normal application behavior must remain unchanged.

Alternative considered: always collecting in memory and only writing on request. Rejected unless implementation can prove zero behavior impact, because unnecessary collection can still change performance or object lifetimes.

### Test Helpers

Add helpers that write and read snapshot JSON in tests. Helpers should keep paths controlled by tests and avoid writing to arbitrary project directories.

Rationale: Snapshot tests need repeatable fixture handling and should not depend on user environment state.

Alternative considered: direct file handling in every test. Rejected because it duplicates schema assumptions and makes updates error-prone.

## Risks / Trade-offs

- Snapshot schema drift -> Mitigation: include schema/version fields and tests that assert required keys.
- Debug capture accidentally changes behavior -> Mitigation: add tests comparing disabled behavior and enabled behavior for the same scenario.
- Snapshot files become too large -> Mitigation: keep metadata optional and restrict required fields to reproducibility data.
- Sensitive or unstable internal data leaks into fixtures -> Mitigation: keep internal metadata optional and avoid requiring timestamps, memory addresses, object reprs, or machine-specific paths.
- Candidate diagnostics are unavailable for some candidates -> Mitigation: require `quality_diagnostic` only when it exists and preserve candidate source information independently.

## Migration Plan

No data migration is required. Snapshot capture is disabled by default and should not affect existing saved documents or normal Clean 2D flows.

Implementation should be introduced behind opt-in gates, then covered by tests before any future algorithmic fixes use the snapshots.

## Open Questions

- What exact environment variable name should be used? The implementation can choose a project-consistent name, but it must be documented and tested.
- Should test fixtures store full snapshot examples or validate generated snapshots structurally? The initial implementation should prefer structural assertions unless a fixture is needed for regression review.
