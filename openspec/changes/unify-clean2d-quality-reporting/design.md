## Context

Clean 2D quality diagnostics are currently scattered across multiple modules (engine.py, safety.py, depiction_candidates.py, Clean2DController). Each module uses its own diagnostic format and vocabulary, causing duplicated logic, fragile tests, and difficulty in understanding overall quality. The document-clean2d-decision-contract change already defined stable states (applied, rejected, preserve-only, no-op, failed-controlled) and reasons (invalid-coordinates, invariant-violation, stereo-risk, boundary-bond-risk, new-crossing-risk, collision-risk, collapsed-ring-risk, excessive-displacement, worse-quality, backend-failure). This change builds on that foundation by defining a stable reporting contract, lightweight helpers, contract tests, and minimal adapters only where existing diagnostic surfaces are already clear.

## Goals / Non-Goals

**Goals:**
- Define a single stable contract shape for Clean 2D diagnostics.
- Map current visual quality metrics to the stable vocabulary from document-clean2d-decision-contract.
- Unify candidate reporting: source module, score/quality value, final state, and rejection reason.
- Preserve internal metadata as diagnostic context while exposing only the stable contract for tests.
- Add minimal adapters/wrappers where existing code already exposes diagnostics clearly.
- Identify conceptual duplication between safety checks, ranking logic, depiction candidate scoring, and controller diagnostics without refactoring it unless the duplication is trivial.

**Non-Goals:**
- Change layout algorithms or RDKit/CoordGen.
- Redesign MolGraph, Canvas, or other core data structures.
- Perform a global refactor of engine.py.
- Change geometry, ranking, or candidate selection behavior.
- Alter user-facing UI messages beyond what is needed to expose the new contract.
- Improve ring handling, aromaticity, macrocycles, or complex imported structures (out of scope for this change).

## Decisions

1. **Lightweight Diagnostic Contract Helpers** - Create a small `clean_2d_quality` module that contains contract helpers and lightweight validation.
   *Rationale:* Keeps contract logic in one place without introducing a broad framework. Prefer dataclass, TypedDict, plain dict helpers, or manual validation. Do not add pydantic or any new dependency unless it already exists in the project. Alternative considered: inline schema definitions per module (rejected due to duplication risk).

2. **Minimal Adapter Surface** - Add wrappers only where an existing diagnostic surface is already clear and can be adapted without changing layout or ranking behavior.
   *Rationale:* Reduces risk and avoids a forced migration across engine.py, safety.py, depiction_candidates.py, and Clean2DController. Alternative considered: migrating all modules to emit only the new contract in this change (rejected due to regression risk).

3. **Internal Metadata Namespace** – Diagnostics will include an `internal` field for module-specific data. Tests can ignore this field when asserting on public contract.
   *Rationale:* Preserves debugging information without polluting the stable contract. Alternative considered: separate internal diagnostic object (rejected due to complexity in passing two objects through call chains).

4. **Reporting-Only Score Normalization** - Quality scores exposed in diagnostics will be normalized to the 0-1 range where 1 is best and 0 indicates failure.
   *Rationale:* Provides stable diagnostic output for tests and reporting. This normalized score MUST NOT change ranking, candidate selection, or layout decisions in this change. Alternative considered: use normalized score for selection immediately (rejected because it changes behavior).

5. **Follow-up for Non-trivial Duplication** - Non-trivial conceptual duplication will be documented as follow-up work instead of removed here.
   *Rationale:* Keeps this change focused on contract, helpers, tests, and minimal adapters. Trivial duplication may be removed if it is local and does not expand scope.

## Risks / Trade-offs

- **Risk:** Existing code may rely on undocumented fields from old diagnostic formats. -> Mitigation: Add adapters only at stable diagnostic boundaries and preserve internal metadata.
- **Risk:** Normalized reporting score could accidentally influence candidate ranking. -> Mitigation: Keep existing ranking and selection inputs unchanged; tests should verify this where feasible.
- **Risk:** Adding wrapper layer may introduce performance overhead in hot paths. -> Mitigation: Keep wrappers minimal and avoid broad module migration.
- **Trade-off:** Internal metadata is preserved but not validated by contract tests, which means bugs in internal data won't be caught automatically. This is acceptable because the stable contract is what matters for consumers and tests.

## Migration Plan

1. Create `clean_2d_quality` module with lightweight contract helpers.
2. Write contract tests validating the stable reporting format before production code changes.
3. Add minimal wrappers only where an existing diagnostic surface is clear.
4. Connect adapters without changing layout, geometry, ranking, or candidate selection.
5. Run targeted and relevant full tests to verify no regressions.
6. Record any non-trivial duplication or broader migration as follow-up work.

## Open Questions

- Should score normalization be configurable per module type (e.g., safety vs depiction)? Currently assuming a global reporting-only 0-1 scale; it must not affect ranking or candidate selection in this change.
