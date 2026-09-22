# Tasks: Clean2D Campaign 4 rigid and multiring layout

## Baseline and audit

- [x] Verify branch, HEAD `346b229`, allowed untracked files, and no unknown changes.
- [x] Read the master campaign policy, architecture catalog, current Clean2D modules, topology contracts, and regression corpus.
- [x] Capture the pre-production baseline in `baseline.md` and `/tmp/chemuson_campaign4_baseline.txt`.

## Rigid-system contract

- [x] Add deterministic JSON-safe rigid-system descriptors using existing multilayer topology.
- [x] Cover fused, spiro, bridged, polycyclic, multiple-rigid-block, and congested attachment metadata.

## Candidate implementation

- [x] Add bounded `rigid_multiring_layout` local orientation candidate.
- [x] Integrate it without changing Campaign 3 global placement or Campaign 5 routing boundaries.
- [x] Apply all existing hard gates and preserve controlled fallback.

## Corpus, evidence, and tests

- [x] Add the twelve required topology-built fixtures or document API limitations.
- [x] Add RED/GREEN tests for descriptors, candidate metadata, determinism, safety, controls, and internal contribution.
- [x] Generate deterministic baseline/current evidence and review every target/non-target family.

## Promotion and lifecycle

- [x] Run focused Campaign 2/3/4 tests, architecture, full suite, compileall, Ruff, diff checks, and strict OpenSpec validation.
- [x] Confirm all promotion gates and record known limitations.
- [x] Commit `Improve Clean2D rigid and multiring layout`.
- [x] Archive this OpenSpec, repair canonical Purpose if needed, validate it, and commit `Archive Clean2D Campaign 4 rigid multiring layout`.
- [x] Do not start Campaign 5.
