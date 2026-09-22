# tighten-clean2d-campaign3-global-placement-safety Specification

## Purpose
Provide explicit, auditable safety gates and deterministic candidate competition for Campaign 3 global block placement.

## Requirements

### Requirement: Global placement has an explicit finite displacement budget

`global_block_placement` SHALL use a finite displacement budget derived from target bond length and topology block counts. A candidate exceeding that budget MUST be rejected with an explicit stable reason and MUST expose both measured displacement and budget in JSON-safe evidence. The exception MUST NOT relax any other hard gate.

#### Scenario: Excessive global displacement is rejected

- **WHEN** a topology-derived global placement exceeds its declared displacement budget
- **THEN** the candidate is rejected, the rejection reason identifies the budget, and the operation retains its safe fallback or another accepted candidate

#### Scenario: Placement within budget can compete

- **WHEN** a global placement is within its declared budget and passes all explicit hard gates
- **THEN** it may compete with other safe candidates and its before/after metrics remain auditable

### Requirement: Global placement evaluates every safety gate explicitly

The global-placement path SHALL explicitly evaluate finite coordinates, selection/invariants, stereo signature, new crossings, collision safety, ring degeneracy, bounding-box sanity, bond-length sanity, and the strategy displacement budget after placement. It MUST NOT infer that later gates passed solely because `is_clean2d_candidate_safe` returned early.

#### Scenario: Early safety return cannot hide a failed gate

- **WHEN** any explicit global-placement hard gate fails after candidate construction
- **THEN** the candidate is rejected with the corresponding reason and its complete hard-gate map is recorded

### Requirement: Complex preserve chooses among safe candidates

The complex-preserve path SHALL compare safe `global_block_placement`, scaffold, and block-unwrap candidates using existing quality metrics and deterministic source tie-breaking. It MUST NOT select global placement solely because that branch appears first, and it MUST preserve the existing preserve-only fallback when no safe candidate remains.

#### Scenario: A safe alternative beats a worse global placement

- **WHEN** global placement and an existing scaffold or unwrap candidate are safe
- **THEN** the deterministic existing quality comparison selects the better candidate and records all safe competing sources

#### Scenario: Protected complex structure remains safe

- **WHEN** global placement exceeds its budget or fails a hard gate for a protected complex graph
- **THEN** no unsafe redraw is selected and the existing preserve-only protection remains observable

### Requirement: Safety evidence covers the Campaign 3 fixtures

The corrective evidence SHALL record max displacement, displacement budget, bounding-box ratio, selected source, competing safe sources, hard-gate results, and rejection reason for all seven Campaign 3 topology fixtures without modifying historical Campaign 3 baseline evidence.

#### Scenario: Seven fixture records remain auditable

- **WHEN** the seven Campaign 3 fixtures are evaluated
- **THEN** each record contains the declared safety fields and remains JSON-safe and deterministic
