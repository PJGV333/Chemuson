# Complete Clean2D Campaign 4 target improvement Specification

## ADDED Requirements

### Requirement: Preserve explicit candidate state semantics

The rigid multiring candidate contract SHALL distinguish construction-time hard-gate safety from post-construction engine acceptance and final selection.

#### Scenario: Boolean hard gates are persisted
- **GIVEN** a topology-derived rigid multiring candidate is constructed
- **WHEN** its metadata is serialized
- **THEN** `hard_gate_checks` is a deterministic JSON-safe mapping whose values are booleans, and `hard_gates_passed` is its conjunction

#### Scenario: Construction does not claim ranking results
- **GIVEN** a candidate before general ranking
- **WHEN** candidate metadata is created
- **THEN** it does not set `accepted_by_engine` or `selected` from hard-gate status; evidence records those states only after engine evaluation

### Requirement: Provide a real target-family improvement

The implementation SHALL produce at least one topology-general rigid/multiring geometry that passes hard gates, measurably improves its declared baseline, survives engine evaluation, and is selected without source-priority bias.

#### Scenario: Promotion target is selected by geometry
- **GIVEN** a target fixture whose baseline needs work
- **WHEN** `run_clean2d_engine` evaluates the candidate set
- **THEN** the selected source is `rigid_multiring_layout`, hard gates pass, target metrics improve, and removing any Campaign 4-specific priority bias does not change the selection

### Requirement: Handle spiro systems by relative ring orientation

The implementation SHALL orient the two ring subspaces of a spiro system independently around their shared center while preserving internal ring geometry, attachments, stereo, and safety gates.

#### Scenario: Spiro sectors reduce the defect
- **GIVEN** a topology-general spiro system with a baseline crossing or rebuild defect
- **WHEN** the rigid multiring candidate is evaluated
- **THEN** the shared spiro center remains fixed, ring sectors become distinct, crossings or quality improve, and all hard gates pass

### Requirement: Correct local exocyclic direction only

For fused or congested systems with substituents, the implementation SHALL change only the first terminal substituent atom when safe, preserve its attachment bond length approximately, and reject transformations that worsen safety or quality.

#### Scenario: Unsafe substitution is preserved
- **GIVEN** a fused system whose outward adjustment worsens a hard metric
- **WHEN** the candidate is evaluated
- **THEN** the candidate is rejected and the safe existing geometry is preserved

### Requirement: Describe true polycyclic and multiple-rigid topology

The descriptor SHALL classify a topology-built system with at least three fused rings as `polycyclic` and SHALL expose `multiple_rigid_blocks` and `rigid_system_count` for linked rigid systems.

#### Scenario: Three fused rings are polycyclic
- **GIVEN** a graph with at least three fused rings
- **WHEN** rigid systems are described
- **THEN** a system reports family `polycyclic` and ring count at least three

#### Scenario: Linked rigid systems are explicit
- **GIVEN** two rigid systems connected by a linker
- **WHEN** rigid systems are described
- **THEN** `multiple_rigid_blocks` is `true` and `rigid_system_count` is at least two

### Requirement: Protect adjacent campaigns and deterministic evidence

The corrective change SHALL preserve Campaign 3 global placement, controls, historical failures, and deterministic JSON evidence.

#### Scenario: Campaign 3 controls remain owned by Campaign 3
- **GIVEN** biphenyl-like, triphenyl-like, branched multiblock, acyclic, aromatic, or stereo-sensitive controls
- **WHEN** the engine runs
- **THEN** no Campaign 4 change reassigns global block placement or introduces a regression

#### Scenario: Evidence separates states
- **GIVEN** baseline and current candidate results
- **WHEN** corrective evidence is written
- **THEN** each target records `hard_gates_passed`, `accepted_by_engine`, `selected`, `source`, metrics before/after, delta, and rejection reason, with `allow_nan=False` JSON serialization and no runtime-dependent canonical field
