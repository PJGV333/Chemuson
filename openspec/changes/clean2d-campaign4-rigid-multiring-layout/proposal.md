# Proposal: Clean2D Campaign 4 rigid and multiring layout

## Why

Campaign 3 provides topology-derived global block placement, but rigid, fused, spiro, bridged, polycyclic, and congested ring systems still lack a general internal-orientation candidate. The existing topology descriptor identifies several ring-system families, yet it does not expose enough attachment geometry to audit or safely improve local rigid orientation.

## What Changes

- Extend the existing multilayer-derived topology evidence with deterministic rigid-system descriptors.
- Add one bounded `rigid_multiring_layout` candidate for local rigid-system orientation and immediate attachment direction.
- Cover fused, spiro, bridged, polycyclic, congested, multiple-rigid-block, and control families with topology-built fixtures.
- Record before/after metrics, hard-gate checks, candidate contribution, and controlled fallback in deterministic JSON evidence.

## Scope

Production scope is limited to the existing Clean2D topology/engine surfaces, the regression fixture catalog, Campaign 4 tests, and this OpenSpec evidence. The implementation reuses `MultilayerChemicalGraph`, `BlockGraph`, motifs, existing geometry and safety metrics, and existing candidate evaluation vocabulary.

## Out of scope

- Campaign 3 global block placement, root selection, parent-child assignment, or global sector allocation.
- Campaign 5 flexible connector routing or long-branch routing.
- Campaign 6 macrocycle/global large-layout search.
- Campaign 7 global candidate ranking or unbounded search.
- Campaign 8 general local polish.
- Molecule-name or fixture-ID production routing.
- GUI, persistence, chemistry identity, architecture catalog, new dependencies, or the four known baseline failures.

## Promotion gates

- Gate A: no new chemical or safety violations.
- Gate B: benzene, simple monocycle, acyclic, and Campaign 3 controls do not regress.
- Gate C: reproducible internal contribution is demonstrated for at least fused/spiro/rigid-multiring target cases where baseline needs work.
- Gate D: the selected improvement is attributable to `rigid_multiring_layout`, not only an external backend.
- Gate E: Campaign 3 global placement remains selected for its promoted controls.
- Gate F: repeated identical runs preserve descriptor, candidate, and result ordering.
- Gate G: unresolved bridged/polycyclic cases use controlled preserve-only fallback.

## Rollback condition

If any hard gate regresses, if target-family improvement is not attributable to the internal candidate, or if controls regress, keep the candidate rejected/experimental and revert the production routing while retaining the evidence and tests.
