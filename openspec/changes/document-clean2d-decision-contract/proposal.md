## Why

Clean 2D already routes through multiple backends, safety gates, ranking heuristics, and canvas integrity checks, but its observable contract is not explicitly specified. Before changing layout algorithms, Chemuson needs a stable contract for what Clean 2D may change, what it must preserve, how modes differ, and how safe failures are reported.

## What Changes

- Define the Clean 2D decision contract as an OpenSpec capability.
- Specify observable inputs, modes, candidate-source vocabulary, result states, rejection reasons, chemical invariants, geometric safety criteria, and non-destructive failure behavior.
- Document that RDKit, CoordGen, `clean2d.v2`, local graph cleaning, scaffold/block strategies, and fallback paths are candidate providers rather than guaranteed outputs.
- Add tasks for contract tests and diagnostic/reporting tests where current support exists.
- Keep this change documentation- and contract-focused; it does not introduce a new layout algorithm or refactor the engine.

Out of scope:

- No new layout algorithm.
- No MolGraph redesign.
- No canvas redesign.
- No public document-format change.
- No RDKit or CoordGen behavior change beyond documenting their role.
- No global refactor of `engine.py`.
- No targeted improvement to rings, macrocycles, imported complex molecules, or other geometry domains yet.

## Capabilities

### New Capabilities

- `clean-2d`: Observable contract for Clean 2D decisions, preservation guarantees, result states, rejection reasons, safety behavior, and documented mode intent.

### Modified Capabilities

- None.

## Impact

- Affected planning area: Clean 2D behavior contract and future regression tests.
- Likely implementation touchpoints for later changes: `src/chemuson/gui/controllers/clean2d_controller.py`, `src/chemuson/clean2d/engine.py`, `src/chemuson/clean2d/safety.py`, `src/chemuson/clean2d/v2.py`, `src/chemuson/clean2d/local_graph_cleaner.py`, `src/chemuson/clean2d/complex_policy.py`, and existing Clean 2D tests.
- No intended changes to runtime behavior in this proposal itself.
