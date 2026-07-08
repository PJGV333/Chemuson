## Why

Clean 2D now has stable decision states, stable reasons, unified `quality_diagnostic`, and opt-in debug snapshots, but complex-structure routing remains unstable. Known failures show that hierarchical and high-risk structures can still be routed through paths that are too local, too global, or insufficiently explained by the current complexity policy.

This change is needed before fixing broader geometry issues so the engine has a testable contract for when complex structures must be preserved, locally repaired, or blocked from global redraw candidates.

## What Changes

- Define a verifiable contract for `classify_clean2d_complexity` and the Clean 2D complex-policy routing decisions.
- Specify which structure signals classify a graph as high risk for global redraw, including hierarchical blocks, multilayer structures, macrocycles, cyclophanes, internal cavities, complex bridges, and imported/complex projections.
- Specify when preserve-only behavior is mandatory for high-risk complex structures.
- Specify when `local_graph_cleaner` is allowed and when it must be rejected for hierarchical or global-risk structures.
- Specify that high-risk complex structures SHALL NOT enter global redraw candidate generation unless the policy explicitly classifies them as safe for that route.
- Require policy decisions to expose stable reasons and optional `quality_diagnostic`/debug snapshot evidence for debugging, without making snapshots part of the algorithm.
- Add tests first for policy classification and routing, especially:
- `test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph`
- `test_complex_engine_does_not_call_global_redraw_candidates`

## Capabilities

### New Capabilities

- `clean-2d-complex-policy`: Defines the contract for Clean 2D complex-structure classification, preserve-only routing, local repair eligibility, global redraw blocking, and diagnostic evidence.

### Modified Capabilities

- None.

## Impact

- Likely affected modules:
- `src/chemuson/clean2d/complex_policy.py`
- `src/chemuson/clean2d/engine.py` only at narrow routing integration points
- `src/chemuson/clean2d/local_graph_cleaner.py` only if eligibility checks need existing signal exposure
- Clean 2D tests covering complex structures, multilayer constraints, and local graph routing

- No expected changes to RDKit, CoordGen, MolGraph, canvas, general candidate ranking, or simple-structure geometry.
- Debug snapshots and `quality_diagnostic` are used as evidence surfaces only; they are not required inputs to policy decisions.
