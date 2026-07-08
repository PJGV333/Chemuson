## Context

Clean 2D currently coordinates UI actions, selection handling, candidate generation, backend-specific layouts, ranking, safety checks, and final canvas application. The key modules are `Clean2DController`, `run_clean2d_engine`, `clean2d.safety`, `clean2d.v2`, `local_graph_cleaner`, `complex_policy`, RDKit isolated/direct helpers, and canvas selection/bond-integrity checks.

The current implementation already contains many protective behaviors, but those behaviors are distributed across modules and are not expressed as a stable observable contract. This change defines the contract first so future algorithmic work can be evaluated against documented behavior instead of implicit heuristics.

## Goals / Non-Goals

**Goals:**

- Define the observable Clean 2D contract for inputs, modes, candidate sources, result states, rejection reasons, invariants, safety gates, and non-destructive behavior.
- Establish stable vocabulary for reporting: `applied`, `rejected`, `preserve-only`, `no-op`, and `failed-controlled`.
- Establish stable rejection reason vocabulary: `invalid-coordinates`, `invariant-violation`, `stereo-risk`, `boundary-bond-risk`, `new-crossing-risk`, `collision-risk`, `collapsed-ring-risk`, `excessive-displacement`, `worse-quality`, and `backend-failure`.
- Identify contract tests and diagnostic/reporting tests to add during implementation.
- Keep future algorithm changes decoupled from this contract work.

**Non-Goals:**

- No new layout algorithm.
- No redesign of `MolGraph`.
- No redesign of canvas or selection internals.
- No public document-format change.
- No RDKit or CoordGen changes beyond documenting their role as candidate providers.
- No global refactor of `engine.py`.
- No targeted improvement to rings, macrocycles, complex structures, or imported depictions in this change.

## Decisions

### Decision: Specify Observable Behavior, Not Backend Internals

The spec defines what Clean 2D must preserve and report, not which backend must win. RDKit, CoordGen, internal templates, `clean2d.v2`, local graph cleanup, scaffold depiction, block unwrap, and fallback strategies remain implementation details exposed only as candidate-source diagnostics.

Alternative considered: specify a preferred backend order as the contract. Rejected because backend availability and quality vary by environment, molecule class, and mode.

### Decision: Use a Stable Result Vocabulary

Clean 2D outcomes will be categorized as `applied`, `rejected`, `preserve-only`, `no-op`, or `failed-controlled`. These categories are intentionally higher-level than existing internal dataclasses so the contract can remain stable while implementation details evolve.

Alternative considered: expose current internal fields directly. Rejected because current classes and messages are implementation-shaped and may change during later refactors.

### Decision: Preserve Chemical Identity Above Visual Improvement

The contract prioritizes invariants over successful geometry changes. Clean 2D must preserve atom IDs, bond IDs, connectivity, bond orders, charges, labels, selection, stereo metadata, and document structure. If no safe candidate exists, the original structure remains intact.

Alternative considered: allow aggressive redraws when visual quality improves. Rejected for this change because professional drawing software must avoid silently changing chemical meaning.

### Decision: Define Minimum Safety Categories Before Numeric Thresholds

The spec names required safety categories such as invalid coordinates, new crossings, collisions, collapsed rings, excessive displacement, worse quality, stereo risk, and boundary-bond risk. Exact numeric thresholds can remain implementation-specific unless a later change formalizes them.

Alternative considered: lock all current thresholds in the spec. Rejected because current values are heuristic and may need calibration during algorithm work.

### Decision: Treat Modes as User Intent

`quick`, `publication`, and `propose` are documented as user-intent modes. `quick` favors safe minimal improvement, `publication` allows more polishing while preserving meaning, and `propose` seeks an alternative safe depiction without requiring the current layout to be bad.

Alternative considered: treat modes only as parameter presets. Rejected because users and tests need a stable behavioral expectation independent of current parameter values.

## Risks / Trade-offs

- Contract may expose gaps in current implementation -> Mitigate by adding tasks that first add tests and diagnostics without requiring algorithm rewrites.
- Stable vocabulary may not map one-to-one to current Spanish/internal messages -> Mitigate with adapter/reporting tests rather than broad message rewrites.
- Over-specification could freeze heuristics too early -> Mitigate by specifying categories and outcomes, not all thresholds.
- Contract tests may initially require small reporting changes -> Mitigate by limiting implementation to reporting/diagnostic surfaces where support already exists.
- Future algorithm work may need new rejection reasons -> Mitigate by allowing later OpenSpec changes to extend the vocabulary explicitly.
