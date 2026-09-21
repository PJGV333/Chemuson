# Design: Clean 2D Campaign Policy

## Overview

The change adds one policy specification and one human roadmap. It is intentionally separate from the existing operational specs: those specs define current behavior; this change defines how future, independently approved campaigns may improve that behavior without weakening the existing contracts.

## Existing contracts reused

The policy uses the existing Clean2D terms and surfaces:

- `run_clean2d_engine`, `Clean2DMode`, `Clean2DCandidate`, `Clean2DResult`, `Clean2DGraphSnapshot`, `assert_clean2d_invariants`, and stable rejection reasons from `src/chemuson/clean2d/engine.py`.
- `Clean2DQualityDiagnostic` and its states, reasons, normalized reporting score, source, and opaque internal metadata from `quality_reporting.py`.
- `Clean2DComplexityProfile`, `preserve_only_required`, `local_repair_allowed`, `global_redraw_allowed`, `internal_route_allowed`, and policy evidence from `complex_policy.py`.
- Opt-in JSON debug snapshots with schema `chemuson.clean2d.debug-snapshot` and version `1`.
- Test-owned corpus, baseline reports, deterministic comparison, metric registry, and review policy under `tests/clean2d_regression/`.

Raw lower-is-better engine scores, higher-is-better normalized reporting scores, and diagnostic geometry metrics are not conflated.

## Master policy model

Every future campaign is evaluated in this order:

1. preserve molecular connectivity and all hard constraints;
2. classify the case by orthogonal size and topology/family attributes;
3. capture a reproducible baseline;
4. generate or evaluate candidates through a declared strategy;
5. reject hard-gate violations before ranking;
6. compare surviving candidates using a metric vector, not a hidden scalar;
7. review target and non-target family effects;
8. promote only through the declared gates, or keep the strategy experimental;
9. preserve evidence and support rollback.

The policy's simple/medium/large/complex-scale posture is normative. It does not choose implementation thresholds for future algorithms.

## Hierarchical direction

The intended future decomposition is connectivity → topology → rigid/semi-rigid blocks → block graph → global block placement → connector orientation → flexible branch routing → internal/local geometry → candidate evaluation → global ranking → local polish. Local polish cannot compensate for an incorrect global/topological decision. This remains architectural direction, not code in this change.

## Hard constraints and metrics

Hard constraints are checked before candidate ranking. The metric vector is diagnostic and may be used to compare candidates only after all hard gates pass. Existing metric definitions and tolerances remain authoritative; new metrics are registered before use and are not assigned arbitrary acceptance thresholds in this change.

## Corpus identity

Case IDs identify fixtures and their meaning; expected quality is separate and may evolve only through an explicit diff. Existing tags such as `baseline`, `known_delicate`, `known_failure`, `complex_policy_guard`, `selection_boundary`, and `stereo_sensitive` remain visible. Topology tags are extensible without renaming existing IDs.

## Observability and determinism

The policy extends existing quality diagnostics, snapshots, baseline reports, and review helpers rather than creating a parallel telemetry system. Normal snapshots remain opt-in. Same input, mode, parameters, and seed must produce reproducible result and candidate ordering within documented tolerances; backend identity, seed, source, and decision evidence are recorded when external behavior is not deterministic.

## Scope guard

The only intended new paths are the active OpenSpec, `docs/clean2d/CAMPAIGN.md`, the contract test, and its baseline record. No production module, architecture catalog, GUI, core, chemio, or existing Clean2D spec is modified because this policy does not add a dependency or alter current behavior. During this phase, that scope claim is verified by manual diff review; no persistent test compares against `HEAD`, `origin/main`, or a merge base.

## Validation

Validate the change and all OpenSpecs in strict mode, run the contract test and the existing Clean2D reporting/corpus/baseline/complex-policy tests, compile the repository, run `git diff --check`, and inspect the final status. Missing local test tools are reported rather than replaced with fabricated results.
