# Proposal: Define Clean 2D Geometry Metric Contracts

## Context
The Clean 2D regression corpus now captures diagnostic geometry metrics for a broad set of molecule families. Those metrics are useful only if their names, units, polarity, availability, and comparison semantics are stable and explicitly documented.

## Objective
Define a formal semantic contract for Clean 2D geometry metrics used by the regression corpus. This change makes metric records reproducible and scientifically interpretable without changing layout, ranking, candidate selection, backends, UI, or document format.

## Scope
- Define stable metric names for metrics already captured by the corpus.
- Document type, unit, polarity, scope, availability, optionality, and diagnostic-only status.
- Define JSON-stable metric record shape for before/after reporting.
- Define explicit floating-point comparison tolerances for tests and serialization checks.
- Add small test-only helpers and tests that validate corpus metric records match the contract.

## Non-Goals
- No geometry improvements.
- No candidate ranking changes.
- No layout algorithm changes.
- No RDKit, CoordGen, or backend changes.
- No UI changes.
- No public document format changes.
- No strict aesthetic gates for current corpus cases.
- No repairs for benzene, fused rings, macrocycles, heteroaromatics, selection boundaries, or multi-block systems.
- No large refactor of `engine.py`.

## Success Criteria
- OpenSpec validates in strict mode before implementation.
- Each captured corpus metric has a stable definition.
- Metric records serialize using JSON-stable primitive values.
- Floating-point comparisons use explicit tolerances.
- Optional and unavailable values are represented consistently.
- Corpus tests continue to pass without failing solely on visually imperfect geometry.
