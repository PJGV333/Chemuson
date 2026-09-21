# Clean 2D Campaign Policy Specification

## Purpose

This specification defines the master policy for future Clean 2D improvement campaigns. It governs experiments, evidence, and promotion without implementing a new layout algorithm or changing an existing Clean2D contract.

## ADDED Requirements

### Requirement: Master policy is documentary and sequential

The campaign policy SHALL govern future Clean2D work as a sequence of independently specified campaigns. This change SHALL NOT implement any campaign or alter production layout behavior.

#### Scenario: Policy change is validated
- **GIVEN** the active campaign-policy change is reviewed
- **WHEN** its files and diff are inspected
- **THEN** it contains policy, roadmap, contract tests, and baseline evidence only
- **AND** no production Clean2D, GUI, core, chemio, persistence, backend, or architecture-catalog behavior is changed.

#### Scenario: Future campaign is proposed
- **GIVEN** a team wants to implement one of Campaigns 1–9
- **WHEN** the work is started
- **THEN** it SHALL use a separate OpenSpec naming the target family, non-target families, expected improvement, hard invariants, observed metrics, baseline input, promotion gates, and rollback condition.

### Requirement: Quality posture is monotonic by complexity

The master policy SHALL apply the following normative posture: a correctly represented simple structure SHALL NOT regress; medium structures SHALL show sustained measurable improvement; large structures SHALL degrade only in a controlled and explainable manner; complex-scale structures SHALL preserve chemistry and use the best safe candidate among the candidates actually evaluated under the declared strategy, search space, and resource budget. This does not claim a mathematical or global optimum. When no safe improvement exists, preserve-only or no-op SHALL be preferred to destructive layout.

For medium structures, “sustained measurable improvement” means reproducible improvement at the target family and campaign level. It does not require every medium structure to improve in every change. An individual medium case that is already good MAY remain unchanged or no-op when it does not regress, the target family shows demonstrable aggregate or distributional improvement, and any individual regression remains visible in review. Campaign 1 SHALL determine which summaries are appropriate; this master policy does not mandate a particular statistical aggregate.

#### Scenario: Simple non-regression
- **GIVEN** a simple corpus case is currently correct
- **WHEN** a future campaign changes Clean2D
- **THEN** the case SHALL not acquire a new contract violation, unsafe geometry, or unexplained result-state regression.

#### Scenario: Medium improvement
- **GIVEN** a target case is medium by size and topology classification
- **WHEN** a campaign claims improvement
- **THEN** the improvement SHALL be reproducible and measurable at the declared target-family and campaign level
- **AND** an already-good individual medium case MAY remain unchanged or no-op without violating this requirement when it does not regress
- **AND** any individual regression SHALL remain visible in review
- **AND** it SHALL not be asserted solely by visual impression.

#### Scenario: Large controlled degradation
- **GIVEN** complexity increases for a large structure
- **WHEN** a candidate is evaluated
- **THEN** chemical identity, safety, and decision observability SHALL remain intact
- **AND** any degradation SHALL be visible in the before/after report and review.

#### Scenario: Complex preserve-only fallback
- **GIVEN** a complex-scale structure has no policy-approved safe improvement
- **WHEN** Clean2D completes
- **THEN** it SHALL preserve the existing geometry or apply only an existing preservation-safe route
- **AND** it SHALL use the existing `preserve-only`, `no-op`, or `failed-controlled` result vocabulary as appropriate.

### Requirement: Hierarchical layout is the architectural direction

Future layout work SHOULD follow the hierarchy molecular connectivity → topological analysis → rigid/semi-rigid block decomposition → block graph → global block placement → connector orientation → flexible branch routing → internal/local geometry → candidate evaluation → global ranking → local polish. Local polish SHALL NOT be treated as a substitute for a wrong topological or global placement decision.

#### Scenario: Local polish is bounded
- **GIVEN** a future campaign adds local polishing
- **WHEN** the polish stage runs
- **THEN** it SHALL preserve the selected global topology and block arrangement unless a separate contract explicitly authorizes otherwise.

### Requirement: Molecule-specific production routing is prohibited

Concrete molecules and regression case IDs SHALL be regression evidence, never production routing rules. Production logic SHALL NOT branch on molecule names, known case IDs, or molecule-specific likeness predicates unless the branch is a demonstrated reusable topological classification.

#### Scenario: Regression case is generalized
- **GIVEN** a case exposes a failure
- **WHEN** a strategy is designed
- **THEN** the strategy SHALL be expressed using reusable signals such as multiblock, macrocycle, fused rings, bridge, congested substitution, flexible connector, or branch competition
- **AND** the case SHALL remain a test fixture.

### Requirement: Corpus taxonomy is orthogonal and extensible

Every corpus case SHALL have a stable case ID and SHALL be describable by a size class and zero or more topology/family tags. Size class SHALL be one of `simple`, `medium`, `large`, or `complex-scale`. Medium is approximately 20–60 graph atoms as an initial objective, but atom count SHALL NOT be the sole complexity criterion.

The extensible tag vocabulary SHALL include, when applicable: `acyclic`, `branched`, `monocycle`, `aromatic`, `multisubstituted-aromatic`, `fused`, `spiro`, `bridged`, `macrocycle`, `multiblock`, `rigid-flexible`, `peptide-like`, `glycoside-like`, `heteroatom-rich`, `charged`, `stereo-sensitive`, `coordination`, `selection-boundary`, `congested`, `known-delicate`, and `known-failure`. Existing corpus tags such as `baseline`, `known_delicate`, `complex_policy_guard`, and `stereo_sensitive` remain compatible.

#### Scenario: Metadata is available
- **GIVEN** a case can be analyzed
- **WHEN** corpus metadata is produced
- **THEN** it SHOULD record `atom_count`, `heavy_atom_count`, `bond_count`, `ring_count`, `connected_components`, `rigid_block_count`, `rotatable_connector_count`, and `macrocycle_count` when available
- **AND** absence of a field in an earlier campaign SHALL NOT invalidate the case ID.

#### Scenario: Case meaning is stable
- **GIVEN** a case ID exists
- **WHEN** its fixture or expected quality changes
- **THEN** the change SHALL be explicit in the diff
- **AND** a known failure SHALL NOT be renamed or silently reclassified to hide a result.

### Requirement: Chemical hard constraints precede aesthetics

A candidate SHALL be rejected before ranking if it violates any existing Clean2D invariant or safety contract. At minimum, the campaign policy treats these as hard constraints: atom IDs, bond IDs, atom count, bond count, element identity, formal charge, bond endpoints, bond order, aromaticity, stereo metadata, applicable selection metadata, finite coordinates, and MolGraph integrity. The existing stable Clean2D result states and rejection reasons SHALL be reused.

#### Scenario: Attractive but unsafe candidate
- **GIVEN** a candidate has a better visual metric
- **WHEN** it changes connectivity, identity, stereo meaning, selection boundary integrity, or finite-coordinate validity
- **THEN** it SHALL be rejected and SHALL NOT win by score.

### Requirement: Hard gates and soft metrics are separate

Hard gates SHALL reject unsafe candidates before ranking. Soft metrics SHALL compare only candidates that pass hard gates. A scalar score MAY summarize a vector but SHALL NOT hide the vector and SHALL NOT compensate for a hard-gate violation.

#### Scenario: Gate ordering
- **GIVEN** several candidates are generated
- **WHEN** ranking is performed
- **THEN** invariant, chemical, stereo, boundary, coordinate, and safety gates SHALL be evaluated first
- **AND** only surviving candidates SHALL be compared by soft metrics.

### Requirement: Quality is represented as a metric vector

Future campaign reports SHALL preserve a vector of diagnostic metrics rather than only a scalar. The vector SHALL include, when applicable: `bond_length_error`, `bond_length_variance`, `bond_angle_penalty`, `atom_collision_count`, `label_collision_count`, `bond_crossing_count`, `ring_distortion`, `ring_degeneracy`, `rigid_block_distortion`, `branch_separation`, `connector_congestion`, `compactness`, `whitespace_balance`, `global_extent`, and `symmetry_preservation`. Reports MAY also include `candidate_count`, `candidate_source`, and `runtime_ms`.

Existing geometry-metric definitions, polarity, optionality, JSON safety, and tolerances remain authoritative. New thresholds SHALL NOT be invented in this policy; Campaign 1 SHALL measure first.

#### Scenario: Metric-only difference
- **GIVEN** two safe candidates differ in diagnostic metrics
- **WHEN** the campaign compares them
- **THEN** the complete vector and metric semantics SHALL remain observable
- **AND** a metric difference SHALL not silently become a chemical gate.

### Requirement: Baselines are captured before algorithm changes

Every future algorithmic campaign SHALL capture a baseline, run the identical corpus, persist the report, make the change, rerun the identical corpus, produce a before/after diff, classify regressions, review relevant cases visually, and explicitly accept or reject the change. Baselines SHALL NOT be updated merely to make tests green; an update SHALL include an explicit justification.

#### Scenario: Before/after review
- **GIVEN** a future strategy changes candidate generation or ranking
- **WHEN** the strategy is evaluated
- **THEN** the report SHALL answer what improved, what regressed, by how much, on which families, why the candidate won, and which hard gates were checked.

### Requirement: Existing observability is extended, not duplicated

Campaign observability SHALL use and extend existing quality diagnostics, debug snapshots, baseline reports, and baseline diff review. Normal debug snapshots SHALL remain opt-in. A diagnostic run SHOULD make strategy, complexity/topology class, candidate sources, candidate count, selected candidate, hard-gate rejections, metrics before/after, decision state, decision reason, and runtime explainable when those values are available.

#### Scenario: Diagnostic evidence is recorded
- **GIVEN** a future campaign evaluates candidates
- **WHEN** a result report is generated
- **THEN** it SHALL preserve existing stable diagnostic states/reasons and source labels
- **AND** policy evidence SHALL be JSON-serializable and observational unless a later OpenSpec changes the contract.

### Requirement: Determinism is measured and reproducible

For the same input, mode, parameters, and seed, Clean2D SHALL aim to reproduce the result and candidate ordering within documented tolerances. If an external backend is nondeterministic, the report SHALL identify the backend, seed when available, candidate source, and final audited decision. Campaign 1 SHALL measure determinism before stricter gates are proposed.

#### Scenario: Deterministic replay
- **GIVEN** a case is replayed with identical inputs and seed
- **WHEN** baseline records are compared
- **THEN** result state, stable reason, candidate source ordering, and selected source SHALL be comparable under the declared tolerance policy.

### Requirement: Performance is measured before SLAs

Campaign 1 SHALL measure wall time, candidate-generation time, and candidate count by size class and topology/family class. The policy SHALL prevent unbounded combinatorial candidate growth. This change SHALL NOT set an arbitrary SLA.

#### Scenario: Candidate search growth
- **GIVEN** a strategy generates multiple candidates
- **WHEN** its performance is measured
- **THEN** candidate count and generation time SHALL be reported
- **AND** uncontrolled growth SHALL block promotion until bounded or explicitly justified.

### Requirement: External backends are candidate sources, not quality definitions

RDKit, CoordGen, or other external backends MAY generate candidates. Backend success SHALL NOT imply layout acceptance. ChemUSON SHALL evaluate every candidate against its own invariants, hard gates, metrics, and review policy.

#### Scenario: Backend output is evaluated by ChemUSON
- **GIVEN** an external backend returns a candidate successfully
- **WHEN** the campaign evaluates that candidate
- **THEN** ChemUSON SHALL apply its own hard gates and metric vector before acceptance
- **AND** backend success alone SHALL not promote the candidate.

### Requirement: Experimental strategies require explicit routing

Experimental strategies MAY run behind explicit routing or benchmark fixtures. They SHALL NOT become the default automatically and SHALL pass the promotion gates before production routing is changed.

#### Scenario: Experimental strategy remains isolated
- **GIVEN** a new strategy has not passed the promotion gates
- **WHEN** it is exercised
- **THEN** it SHALL run only through explicit experimental routing or benchmark fixtures
- **AND** the production default SHALL remain unchanged.

### Requirement: Promotion gates are mandatory

Before a future strategy is promoted, it SHALL pass these gates:

- **Gate A — chemical safety:** zero new invariant or hard-safety violations.
- **Gate B — simple non-regression:** currently correct simple cases do not worsen.
- **Gate C — target-family improvement:** the declared target family improves reproducibly.
- **Gate D — cross-family review:** every non-target regression is identified and justified.
- **Gate E — determinism:** no uncontrolled randomness is introduced.
- **Gate F — performance:** no unjustified performance degradation or unbounded candidate growth.
- **Gate G — manual visual review:** representative before/after cases are reviewed for material changes.

No percentage threshold is fixed here; thresholds must follow evidence and a later OpenSpec.

#### Scenario: Strategy is promoted only after all gates
- **GIVEN** a strategy is proposed for production routing
- **WHEN** promotion is reviewed
- **THEN** Gates A through G SHALL have recorded evidence
- **AND** an unreviewed gate SHALL prevent promotion.

### Requirement: Rollback is safer than molecule exceptions

A strategy SHALL be reverted or kept outside production routing if it violates chemistry, regresses correct simple cases, broadly worsens non-target families, loses determinism, or causes severe unjustified performance degradation. The first response SHALL NOT be a molecule-specific exception. Useful evidence SHALL remain as a test, report, or experiment record.

#### Scenario: Failed strategy is rolled back
- **GIVEN** a promoted strategy introduces a hard violation or broad unexplained regression
- **WHEN** the failure is confirmed
- **THEN** routing SHALL be reverted or the strategy SHALL be kept out of production
- **AND** a molecule-specific exception SHALL not be added as the first response.

### Requirement: Nine campaigns are explicit

The roadmap and future planning SHALL use these ordered campaigns: (1) Benchmark & observability, (2) Topology/decomposition, (3) Medium molecule assembly, (4) Rigid/fused/multiring layout, (5) Flexible connectors and branch routing, (6) Macrocycles and large structures, (7) Global candidate search/ranking, (8) Local polish, and (9) Production acceptance. A campaign SHALL NOT silently combine multiple stages.

#### Scenario: Campaign exit decision
- **GIVEN** a campaign has an exit criterion
- **WHEN** the criterion is evaluated
- **THEN** evidence SHALL be attached to that campaign's OpenSpec
- **AND** the next campaign SHALL not be treated as complete by implication.

### Requirement: Production-ready is evidence-based

Clean2D SHALL be considered production-ready for a campaign only when the declared corpus, hard constraints, deterministic replay, performance evidence, cross-family review, manual visual review, and rollback path are all documented, and the resulting baseline version is identifiable. Visual quality alone SHALL NOT establish production readiness.

#### Scenario: Production readiness requires the evidence bundle
- **GIVEN** a campaign is proposed as production-ready
- **WHEN** the readiness decision is made
- **THEN** corpus, hard-gate, determinism, performance, cross-family, visual-review, baseline, and rollback evidence SHALL be identifiable
- **AND** a visual impression without that evidence SHALL not be sufficient.
