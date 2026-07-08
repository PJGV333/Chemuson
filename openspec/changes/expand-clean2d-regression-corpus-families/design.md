# Design: Expand Clean 2D Regression Corpus Families

## Overview
This change extends the existing test-only corpus rather than changing Clean 2D production behaviour. Builders remain simple, readable, and deterministic. Cases should document current behaviour, including delicate or controlled-failure states, without adding geometry fixes or strict visual acceptance thresholds.

## Case Metadata
The existing case dataclass will be extended with a tuple of tags. Tags are stable strings used by tests and reports:
- `baseline`: ordinary observational case.
- `known_delicate`: current behaviour is useful to record but not yet a geometry requirement.
- `known_failure`: controlled current limitation, if any.
- `complex_policy_guard`: case exercises high-risk or preserve-only policy behaviour.
- `selection_boundary`: case targets a subset and exercises boundary bonds.
- `stereo_sensitive`: case contains stereo metadata that must be preserved.
- `charged`: case contains formal charges that must be preserved.

Existing boolean fields may remain for compatibility with the first corpus change, but tests should treat tags as the primary classification surface.

## New Family Coverage
The expanded registry should keep the original five cases and add a modest set of new cases:
- Substituted aromatics: toluene-like, phenol-like, aniline-like, para-disubstituted benzene-like.
- Additional fused aromatics: naphthalene-like and/or anthracene-like baseline.
- Heteroaromatics: pyridine-like, imidazole-like, furan-like, thiophene-like.
- Charged molecules: ammonium-like, carboxylate-like, pyridinium-like or simple zwitterion-like.
- Stereo-sensitive cases: chiral center with wedge/hash, alkene E/Z-like, and optionally small cyclic stereo.
- Selection-boundary cases: acyclic selection, ring-touching selection, and linker-between-blocks selection.
- Multi-block systems: biphenyl-like, diphenyl-ether-like, triphenyl-like if simple.
- Complex-policy baseline: macrocycle-like or cyclophane/cavity-like if viable with simple builders.

## Assertions
The existing corpus assertions continue to enforce common contracts: valid graph, controlled execution, declared result state, identity preservation, optional metric capture, and JSON-serializable snapshots.

Additional tag-specific assertions will be test-only:
- `charged` cases verify formal charges before and after Clean 2D.
- `stereo_sensitive` cases verify stereo metadata before and after Clean 2D.
- `selection_boundary` cases verify bonds crossing the selected/unselected boundary keep ids, endpoints, order, aromaticity, and stereo metadata.
- `complex_policy_guard` cases verify current high-risk policy behaviour is recorded without requiring an applied global redraw.

## Metrics And Snapshots
Metric capture remains observational and uses existing helpers. Metrics must not become pass/fail criteria for visual quality in this change. Snapshots remain test-owned and JSON-serializable; they do not modify the public document format.

## Gaps And Follow-Ups
If macrocycle, cyclophane, or cavity coverage cannot be represented safely with simple builders, the gap should be documented for a future corpus-expansion change. Strict metric thresholds and visual diffing remain future OpenSpec changes.

## Follow-Up OpenSpec Roadmap
- `define-clean2d-geometry-metric-contracts`: define stable metric names and tolerance semantics after corpus coverage is broader.
- `stabilize-clean2d-deterministic-baselines`: make repeated execution and snapshot comparison deterministic enough for stricter regression gates.
- `improve-clean2d-geometric-quality-v1`: start targeted geometry improvements measured against this corpus.
