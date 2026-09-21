# Tasks: Clean2D Campaign 2 topology/decomposition

## OpenSpec and baseline

- [x] Read the master campaign policy and existing multilayer/decomposition contracts.
- [x] Capture the Campaign 2 baseline before implementation.
- [x] Validate this OpenSpec strictly before implementation.

## Deterministic decomposition evidence

- [x] Add a JSON-safe topology decomposition summary using existing multilayer and block graph objects.
- [x] Preserve stable block IDs, block kinds, atom membership, anchors, motif IDs, and metadata.
- [x] Expose connector/edge evidence without adding a parallel molecular graph.
- [x] Canonicalize connected components, counts, and all collections deterministically.
- [x] Keep geometry, candidate generation, ranking, local repair, and backend routing unchanged.

## Tests and gates

- [x] Add contract tests for medium/large block decomposition and connector evidence.
- [x] Add determinism and JSON-serialization tests.
- [x] Re-run focused topology tests, architecture tests, and the full suite.
- [x] Run compileall, Ruff, strict OpenSpec validation, and diff checks.
- [x] Verify no unauthorized production or architecture-catalog paths changed.

## Promotion

- [x] Commit exactly the Campaign 2 scope as `Build deterministic Clean2D topology decomposition evidence`.
- [x] Archive this OpenSpec only after all gates pass and the commit is verified.
