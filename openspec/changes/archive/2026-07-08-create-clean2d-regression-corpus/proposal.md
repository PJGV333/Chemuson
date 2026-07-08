# Proposal: Create Clean 2D Regression Corpus

## Context
Clean 2D needs a stable scientific baseline before targeted geometry improvements are attempted. Existing tests cover individual behaviours, but they do not provide a named corpus that can be reused to measure regressions across representative molecule families and edge cases.

## Objective
Create regression infrastructure for Clean 2D: a stable registry of named molecule cases, contract assertions, JSON-serializable debug snapshots, and optional metric capture. This change records current behaviour and protects chemical identity invariants; it does not improve geometric output.

## Scope
- Add a regression case registry with stable names and metadata: family, mode, target, expected contract states, and known-delicate markers.
- Add molecule builders or fixtures for a small initial set of representative cases.
- Add contract assertion helpers for uncontrolled exceptions, chemical identity preservation, expected result states, and debug snapshot serializability.
- Capture optional metrics only through existing helpers where available, or through minimal test-only adapters when a small helper is missing.
- Treat `quality_reporting.py` as a stable reporting, normalization, and diagnostic layer, not as a geometry metrics engine.

## Non-Goals
- No new layout algorithms.
- No large refactor of `engine.py`.
- No UI changes.
- No direct improvements to tetrandrine, naphthalene, macrocycles, or smart propose behaviour.
- No RDKit or CoordGen changes.
- No public document format changes.
- No candidate ranking or geometry generation changes.

## Success Criteria
- OpenSpec validates in strict mode.
- The regression corpus runs through the Clean 2D engine without uncontrolled exceptions.
- Each case declares family, mode, target, and expected contract states.
- Each case preserves chemical identity invariants.
- Each case can produce a JSON-serializable debug snapshot.
- Known-delicate cases can record current behaviour without forcing geometry fixes.
