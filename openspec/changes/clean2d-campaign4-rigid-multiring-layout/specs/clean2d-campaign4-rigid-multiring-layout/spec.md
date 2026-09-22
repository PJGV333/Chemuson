# Specification: Clean2D Campaign 4 rigid and multiring layout

## ADDED Requirements

### Requirement: Rigid systems are described through existing topology

Campaign 4 SHALL derive a deterministic, JSON-safe descriptor from the existing `MultilayerChemicalGraph`/`BlockGraph` and molecular graph. Each rigid-system record SHALL expose member atom IDs, member bond IDs, ring membership, shared atoms and bonds, external attachment atoms and neighbors, centroid, principal orientation when derivable, attachment vectors, external substituent count, and local congestion. The descriptor SHALL use topology and geometry properties, never molecule names or fixture IDs.

#### Scenario: Descriptor is stable and auditable

- **GIVEN** the same graph, selected atoms, and coordinates
- **WHEN** the rigid-system descriptor is generated twice
- **THEN** the records, ordering, IDs, and JSON serialization are identical
- **AND** all unavailable geometric values are explicit `null` rather than NaN or infinity.

#### Scenario: Families are classified generally

- **GIVEN** ring systems with one ring, fused overlap, one-atom spiro overlap, bridge paths, or multiple connected rigid systems
- **WHEN** topology is described
- **THEN** the descriptor exposes the applicable `monocycle`, `fused`, `spiro`, `bridged`, `polycyclic`, or `multiple_rigid_blocks` family without identity-specific routing.

### Requirement: Campaign 4 emits one bounded internal rigid-layout candidate

The Clean2D engine SHALL be able to emit a candidate source named `rigid_multiring_layout` for reusable rigid-system signals. Its strategy MAY rotate a rigid system around its existing centroid and choose immediate outward attachment directions, but SHALL NOT translate complete blocks, change Campaign 3 root/parent assignments or global sector allocation, route flexible branches, or perform unbounded global search.

#### Scenario: Local rigid orientation is attempted

- **GIVEN** a fused, spiro, bridged, polycyclic, congested, or multiple-rigid-block graph requiring local work
- **WHEN** Campaign 4 candidate generation runs
- **THEN** it emits at most the bounded candidate count declared by the strategy
- **AND** metadata records strategy, rigid-system count/types, affected atoms, attachments, and before/after quality.

#### Scenario: Simple and already-good controls are not forced

- **GIVEN** regular benzene, a simple monocycle, an acyclic graph, or a rigid case already classified as good
- **WHEN** candidate generation runs
- **THEN** Campaign 4 does not force a geometry change and existing control behavior remains available.

### Requirement: Rigid-layout candidates pass existing hard gates

A Campaign 4 candidate SHALL preserve atom and bond identity, bond endpoints/order/aromaticity, element and charge identity, stereo signature, finite coordinates, selection integrity, no new crossings, collision safety, ring degeneracy, bond-length sanity, bounding-box sanity, and a finite local displacement budget. Hard gates SHALL be evaluated before quality comparison; no gate may be relaxed for Campaign 4.

#### Scenario: Unsafe rigid orientation is rejected

- **GIVEN** a proposed local orientation violates any hard gate
- **WHEN** the candidate is evaluated
- **THEN** it is rejected with a stable reason and a complete JSON-safe hard-gate map
- **AND** the engine keeps an accepted existing candidate or controlled preserve-only fallback.

### Requirement: Safe local improvement is attributable to the internal candidate

A `rigid_multiring_layout` candidate SHALL be accepted only when it passes hard gates and improves a declared relevant ring/attachment/visual metric relative to the selected baseline, unless the existing policy explicitly preserves a good no-op. Evidence SHALL distinguish candidate emission, candidate acceptance, final selection, and external-backend fallback.

#### Scenario: Internal contribution is observable

- **GIVEN** a target-family case where the baseline needs work
- **WHEN** Campaign 4 is evaluated
- **THEN** the report identifies whether `rigid_multiring_layout` was emitted, accepted, and selected
- **AND** its before/after metrics and delta are available independently of any external backend.

### Requirement: Campaign 3 and Campaign 5 boundaries remain intact

Campaign 4 SHALL consume Campaign 3 block-placement metadata when present but SHALL NOT change global block placement, root block, parent-child assignment, global sector allocation, flexible connector torsion, or long-branch routing. Unresolved bridged/polycyclic cases MAY remain preserve-only.

#### Scenario: Cross-campaign boundaries are preserved

- **GIVEN** a Campaign 3 medium multiblock case or a graph with a long flexible branch
- **WHEN** Campaign 4 candidate generation runs
- **THEN** Campaign 3 global placement and branch-routing behavior remain unchanged
- **AND** any unresolved rigid case uses the existing safe fallback.

### Requirement: Promotion evidence covers target and control matrix

Campaign 4 SHALL provide deterministic baseline/current evidence for regular benzene, simple monocycle, fused aromatic and non-aromatic systems where supported, spiro, bridged, larger polycyclic, fused substitutions, spiro substitution, two rigid blocks with linker, and congested attachments. Evidence SHALL omit `runtime_ms` from deterministic comparisons and SHALL report limitations instead of inventing unsupported chemistry.

#### Scenario: Matrix and gates are reviewable

- **GIVEN** the required fixture matrix is evaluated twice
- **WHEN** evidence is compared
- **THEN** candidate ordering and stable outcomes are identical
- **AND** chemical safety, controls, internal contribution, Campaign 3 preservation, determinism, and controlled fallback are separately reviewable.
