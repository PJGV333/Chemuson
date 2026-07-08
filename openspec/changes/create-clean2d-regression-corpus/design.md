# Design: Clean 2D Regression Corpus

## Overview
The corpus is a test-only baseline for measuring Clean 2D stability. It records representative molecule inputs, the expected Clean 2D contract state, identity invariants, optional diagnostics, and serializable debug snapshots. It intentionally avoids changing geometry generation, ranking, UI behaviour, RDKit/CoordGen integration, or document format.

## A. Regression Case Registry
The registry defines stable named cases using a frozen dataclass. Each case declares:
- `name`: stable pytest id and debug snapshot name.
- `family`: grouping such as `acyclic`, `aromatic`, `fused`, `partial-selection`, or `dirty-import`.
- `builder`: test-only molecule builder or fixture loader.
- `mode`: Clean 2D mode used by the engine.
- `target`: whole molecule or an explicit selection target.
- `expected_states`: accepted Clean 2D contract states for current behaviour.
- `known_delicate` and optional `known_failure`: markers for current behaviour that should be recorded without forcing a geometry fix.

The registry is intentionally data-oriented. It should not contain geometry algorithms.

## B. Molecule Builders / Fixtures
Builders create valid molecule graphs with stable atom ids, bond ids, element symbols, bond orders, and coordinates sufficient for Clean 2D execution. Initial cases should be small and representative rather than exhaustive.

Fixtures may be used where existing project files already encode useful inputs. New public document format fields must not be introduced.

## C. Contract Assertions
Contract helpers validate that each case:
- Builds a valid molecule graph.
- Runs through Clean 2D without uncontrolled exceptions.
- Produces an accepted contract state.
- Preserves chemical identity invariants such as atom ids, bond ids, element symbols, bond endpoints, bond orders, and stable stereo metadata where present.
- Produces a JSON-serializable debug snapshot.

The assertion layer should stay focused on test contracts. It must not become a replacement Clean 2D engine.

## D. Optional Metric Collection
Metric collection is optional and observational. Metrics such as crossings, collisions, ring degeneracy, displacement, or visual quality may be captured only by reusing existing helpers where they already exist.

If a small helper is missing, this change may add a minimal test-only adapter. `quality_reporting.py` remains a stable reporting, normalization, and diagnostic layer; it should not become a heavy geometry metrics module.

Metrics captured by this corpus must not influence candidate ranking, geometry generation, layout decisions, or expected pass/fail state unless a future OpenSpec change explicitly defines that contract.

## E. Optional Debug Snapshot Generation
Each case should be able to emit a JSON-serializable debug snapshot containing case metadata, pre/post identity signatures, result state, optional reason, optional diagnostics, and optional metrics. Snapshot generation is for diagnostics and future visual comparisons only.

The snapshot schema is test-owned and must not alter the public document format.

## F. Known-Failure / Known-Delicate Case Handling
Known-delicate cases record current behaviour for scientifically useful inputs that may have imperfect geometry today. They may allow multiple current contract states when the behaviour is intentionally being observed rather than fixed.

Known-failure cases, if added, must distinguish controlled failures from uncontrolled exceptions. They should document the current limitation and point to a follow-up OpenSpec roadmap item rather than hiding real crashes.

## Initial Coverage
The first implementation should include a minimal set of stable cases:
- Simple acyclic molecule.
- Simple aromatic ring.
- Fused aromatic or fused-ring representative, marked delicate if needed.
- Partial-selection or boundary-sensitive case.
- Dirty-coordinate valid graph case.

Complex tetrandrine, macrocycle, naphthalene quality improvements, and smart-propose improvements remain future work unless represented only as known-delicate observational cases.

## Follow-Up OpenSpec Roadmap
- `expand-clean2d-regression-corpus-families`: add broader molecule families and fixture coverage.
- `define-clean2d-geometry-metric-contracts`: specify stable metrics and thresholds after helpers are mature.
- `improve-clean2d-geometric-quality-v1`: targeted geometry improvements measured against this corpus.
- `add-clean2d-visual-regression-diffs`: optional visual diff tooling built on snapshots.
