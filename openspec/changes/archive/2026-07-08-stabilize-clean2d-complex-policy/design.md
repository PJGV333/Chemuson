## Context

Clean 2D already exposes stable decision states, stable reasons, quality diagnostics, and opt-in debug snapshots. Those surfaces make decisions observable, but the policy that routes complex structures between preserve-only, local repair, internal candidates, and global redraw candidates is still not explicit enough.

The current risk is algorithmic routing rather than serialization or reporting: hierarchical blocks, multilayer structures, macrocycles, cyclophanes, internal cavities, and imported complex projections need a stable policy before broader layout fixes are attempted. The implementation should remain small and testable, focused on `complex_policy.py` and narrow routing checks in `engine.py`.

## Goals / Non-Goals

**Goals:**

- Define a verifiable contract for `classify_clean2d_complexity` and complex-policy outputs.
- Classify high-risk complex structures deterministically from existing structural, multilayer, and quality signals.
- Make preserve-only mandatory for structures classified as too risky for redraw.
- Gate `local_graph_cleaner` so it is used only for structures where local repair is policy-approved.
- Prevent high-risk complex structures from reaching global redraw candidate generation unless the policy explicitly allows that route.
- Preserve existing safety gates and stable reasons.
- Use `quality_diagnostic` and opt-in debug snapshots as evidence for policy decisions.
- Add tests first for the two in-scope failures before implementation changes.

**Non-Goals:**

- Do not fix naphthalene or fused aromatic systems in this change.
- Do not fix propose mode in this change.
- Do not implement a new layout engine.
- Do not perform a global `engine.py` refactor.
- Do not change RDKit, CoordGen, MolGraph, or canvas.
- Do not change geometry for simple structures.
- Do not weaken existing Clean 2D safety gates.
- Do not use debug snapshots as required algorithm inputs.

## Decisions

### Policy-First Classification

`classify_clean2d_complexity` should return or expose enough stable policy information to explain whether a structure is eligible for preserve-only, local repair, internal candidates, or global redraw candidates. The policy should use existing graph, multilayer, and quality signals instead of recomputing layout candidates.

Rationale: routing decisions become testable without changing candidate generation or ranking.

Alternative considered: patching individual failing tests in `engine.py`. Rejected because it would hide policy drift and create more special cases.

### High-Risk Signals Block Global Redraw

Hierarchical blocks, macrocycles, cyclophanes, internal cavities, complex bridged systems, and imported complex projections should be treated as high-risk unless the policy has an explicit safe route. High-risk classification should block global redraw candidate generation.

Rationale: global redraw candidates can collapse or reinterpret complex topology. Blocking them at the policy boundary avoids accidental calls to RDKit/CoordGen-style redraw paths for structures where preservation is safer.

Alternative considered: generate all candidates and reject unsafe ones later. Rejected for high-risk structures because the route itself is disallowed by the intended policy and tests must be able to assert candidate generators were not called.

### Local Repair Eligibility Is Explicit

`local_graph_cleaner` should be allowed only when the policy classifies the structure as suitable for local repair. Hierarchical block structures and structures with global-risk signals should not be selected into local repair unless a specific local defect policy allows it.

Rationale: local repair can be correct for some complex structures, but it must not become a default substitute for global redraw blocking.

Alternative considered: always try local repair before preserve-only. Rejected because hierarchical cases can select the wrong local graph route.

### Preserve-Only Is a Stable Outcome

For high-risk structures without an approved repair route, Clean 2D should preserve the current geometry with stable state/reason evidence rather than attempting global redraw or unstable local repair.

Rationale: preserve-only is safer than geometry collapse for complex imported structures and aligns with the existing decision contract.

Alternative considered: return failed-controlled for all high-risk structures. Rejected because some high-risk cases are valid drawings where preserving the projection is the intended safe outcome.

### Diagnostics Are Evidence, Not Inputs

`quality_diagnostic` and opt-in debug snapshots should record policy reasons, candidate-route decisions, and final state where available, but policy decisions must not require snapshots to be enabled.

Rationale: debugging must be reproducible while normal behavior remains independent of snapshot capture.

Alternative considered: making snapshot contents part of policy tests. Rejected because snapshots are opt-in diagnostic artifacts, not algorithmic dependencies.

## Risks / Trade-offs

- Risk: Over-blocking complex structures prevents legitimate repairs -> Mitigation: tests should cover preserve-only and local-repair eligibility separately, and simple structures must remain unchanged.
- Risk: Policy data shape grows too broad -> Mitigation: expose only stable route flags/reasons needed for tests and diagnostics.
- Risk: Routing changes accidentally affect candidate ranking -> Mitigation: do not change rank ordering; only block or allow routes according to policy contract.
- Risk: Existing safety gates are bypassed -> Mitigation: preserve safety validation after any allowed candidate path.
- Risk: Fix requires broad `engine.py` redesign -> Mitigation: stop and record a follow-up rather than refactor globally.

## Migration Plan

No saved-data migration is required. The change affects Clean 2D runtime routing only.

Implementation should start with failing/pending tests for policy classification and route blocking, then add the smallest policy and routing changes necessary to satisfy them.

## Open Questions

- What stable reason labels should distinguish preserve-only due to high global-redraw risk versus local-repair rejection? The implementation should reuse existing stable vocabulary where possible and add only minimal policy-specific metadata if needed.
