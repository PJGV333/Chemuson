# clean-2d Specification

## Purpose
TBD - created by archiving change document-clean2d-decision-contract. Update Purpose after archive.
## Requirements
### Requirement: Clean 2D accepts explicit operation inputs
Clean 2D SHALL define its operation input as a molecular graph, a target atom set or whole-structure target, the current coordinates, the active mode, and the requested target bond length.

#### Scenario: Whole-structure clean input
- **GIVEN** a document with atoms and structural bonds and no active structural selection
- **WHEN** Clean 2D is invoked
- **THEN** the operation target SHALL be the complete molecular structure

#### Scenario: Selection clean input
- **GIVEN** a document with a selected connected or disconnected atom set
- **WHEN** Clean 2D is invoked
- **THEN** the operation target SHALL be the selected atoms and their structural bonds, with boundary bonds considered for integrity validation

### Requirement: Clean 2D modes have documented intent
Clean 2D SHALL expose documented mode intent for `quick`, `publication`, and `propose`.

#### Scenario: Quick mode intent
- **GIVEN** a user invokes quick Clean 2D
- **WHEN** candidate selection is performed
- **THEN** the operation SHALL prefer safe minimal cleanup over aggressive redrawing

#### Scenario: Publication mode intent
- **GIVEN** a user invokes publication Clean 2D
- **WHEN** candidate selection is performed
- **THEN** the operation SHALL allow stronger geometric polishing than quick mode while preserving chemical identity

#### Scenario: Propose mode intent
- **GIVEN** a user invokes propose Clean 2D
- **WHEN** candidate selection is performed
- **THEN** the operation SHALL seek a safe alternative depiction and SHALL NOT require that the current depiction is chemically invalid

### Requirement: Clean 2D candidate sources are reportable
Clean 2D SHALL treat candidate sources as reportable origins of proposed coordinates without guaranteeing that any specific source wins.

#### Scenario: Candidate source recorded
- **GIVEN** a Clean 2D run evaluates one or more candidate coordinate sets
- **WHEN** a result is produced
- **THEN** each reportable candidate SHALL have a source label suitable for diagnostics

#### Scenario: Backend failure is controlled
- **GIVEN** a backend such as RDKit, CoordGen, local graph cleaning, scaffold depiction, block unwrap, internal templates, `clean2d.v2`, or fallback polishing cannot produce a usable candidate
- **WHEN** other candidates or failure handling continue
- **THEN** the failed backend SHALL be represented as `backend-failure` or equivalent diagnostic detail without mutating the original structure

### Requirement: Clean 2D result states are stable
Clean 2D SHALL classify observable outcomes using the stable result states `applied`, `rejected`, `preserve-only`, `no-op`, and `failed-controlled`.

#### Scenario: Applied result
- **GIVEN** a safe candidate improves or validly changes the targeted coordinates
- **WHEN** Clean 2D commits the candidate
- **THEN** the result state SHALL be `applied`

#### Scenario: Rejected result
- **GIVEN** a candidate violates invariants or safety criteria
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with a stable rejection reason

#### Scenario: Preserve-only result
- **GIVEN** a structure is classified as too complex or risky for redraw but eligible for preservation-oriented handling
- **WHEN** Clean 2D completes without a destructive redraw
- **THEN** the result state SHALL be `preserve-only` if the operation intentionally preserved the original layout or applied only preservation-safe repair

#### Scenario: No-op result
- **GIVEN** the current depiction is already acceptable or no meaningful coordinate movement is needed
- **WHEN** Clean 2D completes without coordinate changes
- **THEN** the result state SHALL be `no-op`

#### Scenario: Failed-controlled result
- **GIVEN** no safe candidate can be produced or all candidates fail controlled checks
- **WHEN** Clean 2D completes
- **THEN** the result state SHALL be `failed-controlled` and the original structure SHALL remain intact

### Requirement: Clean 2D preserves chemical identity
Clean 2D SHALL preserve chemical identity and document structure while changing only accepted target coordinates.

#### Scenario: Atom and bond identity preserved
- **GIVEN** a molecule before Clean 2D
- **WHEN** Clean 2D returns any result state
- **THEN** atom IDs, bond IDs, atom count, bond count, and connectivity SHALL be preserved

#### Scenario: Chemical attributes preserved
- **GIVEN** atoms and bonds with orders, aromatic flags, charges, labels, isotopes, radicals, query flags, mapping, and display metadata
- **WHEN** Clean 2D returns any result state
- **THEN** those non-coordinate chemical and display attributes SHALL be preserved

#### Scenario: Stereo metadata preserved
- **GIVEN** atoms or bonds with stereochemical metadata or wedge/hash styles
- **WHEN** Clean 2D returns any result state
- **THEN** the stereochemical metadata and bond styles SHALL be preserved

#### Scenario: Selection preserved
- **GIVEN** a user has an active atom or bond selection before Clean 2D
- **WHEN** Clean 2D completes
- **THEN** the selection SHALL be preserved unless a later explicit user action changes it

### Requirement: Clean 2D is non-destructive on unsafe output
Clean 2D SHALL leave the original structure intact when no safe candidate exists.

#### Scenario: No safe candidate
- **GIVEN** all candidate coordinate sets are unsafe, invalid, or unavailable
- **WHEN** Clean 2D completes
- **THEN** no atom coordinates or non-coordinate molecular data SHALL be changed

#### Scenario: Candidate rejected after evaluation
- **GIVEN** a candidate initially appears available
- **WHEN** later invariant, safety, stereo, or boundary validation rejects it
- **THEN** the original structure SHALL remain unchanged

### Requirement: Clean 2D enforces minimum geometric safety categories
Clean 2D SHALL reject candidates that violate minimum geometric safety categories.

#### Scenario: Invalid coordinates rejected
- **GIVEN** a candidate has missing, non-finite, or non-numeric coordinates for targeted atoms
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `invalid-coordinates`

#### Scenario: New crossing risk rejected
- **GIVEN** a candidate introduces new structural bond crossings compared with the original depiction
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `new-crossing-risk`

#### Scenario: Collision risk rejected
- **GIVEN** a candidate creates unacceptable non-bonded atom collisions or atom-bond collisions
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `collision-risk`

#### Scenario: Collapsed ring risk rejected
- **GIVEN** a candidate collapses a ring or makes ring geometry degenerate
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `collapsed-ring-risk`

#### Scenario: Excessive displacement rejected
- **GIVEN** a candidate moves targeted atoms beyond the mode-appropriate safe displacement envelope
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `excessive-displacement`

#### Scenario: Worse quality rejected
- **GIVEN** a candidate does not improve a depiction that requires cleanup or makes accepted quality metrics worse
- **WHEN** Clean 2D evaluates the candidate
- **THEN** the candidate SHALL be rejected with `worse-quality`

### Requirement: Clean 2D enforces chemical and canvas safety categories
Clean 2D SHALL reject candidates that risk chemical identity, stereochemical meaning, or canvas boundary integrity.

#### Scenario: Invariant violation rejected
- **GIVEN** a candidate changes anything other than allowed target coordinates
- **WHEN** Clean 2D evaluates invariants
- **THEN** the candidate SHALL be rejected with `invariant-violation`

#### Scenario: Stereo risk rejected
- **GIVEN** a candidate risks changing the visual or metadata interpretation of stereochemistry
- **WHEN** Clean 2D evaluates stereo preservation
- **THEN** the candidate SHALL be rejected with `stereo-risk`

#### Scenario: Boundary bond risk rejected
- **GIVEN** Clean 2D targets a selection with bonds crossing the selection boundary
- **WHEN** a candidate would distort or disconnect those boundary bonds beyond accepted limits
- **THEN** the candidate SHALL be rejected with `boundary-bond-risk`

### Requirement: Clean 2D reports rejection reasons using stable vocabulary
Clean 2D SHALL map controlled candidate rejection and backend failure cases to stable reason vocabulary for tests and diagnostics.

#### Scenario: Stable rejection reason emitted
- **GIVEN** a candidate is rejected for invalid coordinates, invariant violation, stereo risk, boundary bond risk, new crossing risk, collision risk, collapsed ring risk, excessive displacement, worse quality, or backend failure
- **WHEN** Clean 2D reports the candidate outcome
- **THEN** the report SHALL include the corresponding stable reason value

#### Scenario: Implementation detail may be included
- **GIVEN** a stable rejection reason is reported
- **WHEN** additional backend-specific or heuristic detail is available
- **THEN** Clean 2D MAY include diagnostic detail while preserving the stable reason value

### Requirement: Clean 2D application is undoable when coordinates change
Clean 2D SHALL apply accepted coordinate changes as a single undoable operation.

#### Scenario: Applied coordinate change can be undone
- **GIVEN** Clean 2D applies a candidate that changes coordinates
- **WHEN** the user invokes undo once
- **THEN** the target coordinates SHALL return to their pre-clean positions

#### Scenario: Rejected candidate does not affect undo stack
- **GIVEN** Clean 2D rejects all candidates
- **WHEN** the operation completes
- **THEN** no coordinate-changing undo command SHALL be pushed for rejected candidates

### Requirement: Clean 2D improves simple aromatic layout selection
Clean 2D SHALL be able to select a safe, measurably better layout candidate for simple isolated or substituted aromatic molecules without weakening chemical identity or complex-structure protections.

#### Scenario: Simple aromatic candidate improves measurably
- **GIVEN** a simple isolated or substituted aromatic molecule with an available safe candidate that improves at least one tracked geometry or visual metric without worsening crossing or ring-degeneracy safety
- **WHEN** Clean 2D selects and applies a candidate
- **THEN** at least one tracked simple-aromatic corpus case SHALL show a measurable improvement such as higher visual score, lower bond-length error, lower angle deviation, or safer nonbonded spacing

#### Scenario: Chemical identity remains preserved
- **GIVEN** Clean 2D improves a simple aromatic layout
- **WHEN** the result is compared with the input molecule
- **THEN** atom IDs, bond IDs, atom count, bond count, elements, charges, bond endpoints, bond orders, aromaticity, stereo metadata, and selection metadata SHALL remain unchanged

#### Scenario: Protected complex behavior remains unchanged
- **GIVEN** a high-risk complex-policy molecule or protected tetrandrine-like molecule
- **WHEN** Clean 2D runs after the simple aromatic improvement
- **THEN** the result state and protection behavior SHALL NOT change because of the simple aromatic improvement

#### Scenario: Baseline report changes remain reviewable
- **GIVEN** before and after Clean 2D baseline reports for the simple aromatic improvement
- **WHEN** the reports are compared and reviewed
- **THEN** changed cases SHALL be attributable to expected simple-aromatic behavior or observational geometry metrics
- **AND** the change SHALL NOT introduce unexpected contract-risk outside the intended family

