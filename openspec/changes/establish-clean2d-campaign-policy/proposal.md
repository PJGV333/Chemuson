# Proposal: Establish the Clean 2D Campaign Policy

## Why

Clean 2D is a strategic ChemUSON capability, but the next improvements must be made as a controlled sequence rather than as isolated molecule-specific repairs. Existing contracts already define chemical preservation, stable result states, diagnostic reporting, complex-policy routing, snapshots, regression cases, metrics, baselines, and diff review. A master policy is needed to connect those contracts and govern the future campaign.

## What Changes

- Establish a normative master policy for nine sequential Clean 2D campaigns.
- Define the simple/medium/large/complex-scale quality posture.
- Define orthogonal corpus taxonomy and stable case identity rules.
- Define hard constraints, hard gates, diagnostic metric vectors, baseline comparison, determinism, observability, performance, visual review, promotion, and rollback contracts.
- Define hierarchical layout as architectural direction without implementing algorithms.
- Require a separate OpenSpec for every future campaign.
- Add a human roadmap at `docs/clean2d/CAMPAIGN.md` and a structural contract test.

## Scope

This is a documentation and contract change. It adds no production Clean2D algorithm, routing rule, threshold, backend, GUI behavior, persistence behavior, MolGraph behavior, or new production module.

## Non-goals

- No Campaign 1 implementation.
- No geometry, topology, ranking, candidate-generation, safety, or complex-policy changes.
- No molecule-specific production routing.
- No new metric thresholds before sufficient baseline evidence exists.
- No changes to existing Clean2D, regression, baseline, snapshot, or geometry-metric contracts.
- No OpenSpec archive or merge.

The “no production changes” assertion for this phase is performed through the reviewed diff and commit scope. It is not encoded as a persistent test against `HEAD`, `origin/main`, or a merge base, because this OpenSpec remains active while later campaigns add their own changes.

## Compatibility

The policy reuses the existing vocabulary and contracts from `clean-2d`, `clean-2d-quality-reporting`, `clean-2d-complex-policy`, `clean-2d-debug-snapshots`, `clean-2d-regression-corpus`, `clean-2d-baseline-reports`, `clean-2d-baseline-diff-review`, `clean-2d-geometry-metrics`, and deterministic-baseline changes. Existing diagnostic scores and geometry metrics remain observational unless a later change explicitly promotes them.

## Success criteria

- The active OpenSpec validates in strict mode.
- The roadmap contains all nine campaigns and the required policy sections.
- The master spec distinguishes hard gates from soft metrics and preserves the existing stable vocabularies.
- A structural test verifies scope guards and the main normative policy clauses.
- No production Clean2D file is changed.
