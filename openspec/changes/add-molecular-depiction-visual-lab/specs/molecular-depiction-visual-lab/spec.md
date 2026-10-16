# Molecular Depiction Visual Lab Specification

## Purpose
Provide an isolated Qt development tool for visual calibration of ChemUSON molecular depiction using the production renderer without modifying chemical data or default production behavior.

## ADDED Requirements

### Requirement: Visual lab renders through the production molecular renderer
The lab SHALL display production `ChemusonCanvas`, `AtomItem`, and `BondItem` render output in its Current and Working panels. It SHALL NOT reproduce bonds as independent decorative lines or call Clean2D automatically. Current and Working SHALL represent independent copies of the same input graph and coordinate layout.

#### Scenario: Same structure in both panels
- **GIVEN** a stable gallery case and a Working style profile
- **WHEN** Current and Working are rendered
- **THEN** both canvases contain equivalent atom/bond identity, coordinates, connectivity, orders, aromaticity, charge and stereo metadata
- **AND** only the Working rendering style differs.

### Requirement: Current profile preserves existing defaults
The Current profile SHALL be read-only and SHALL reproduce the exact `CHEMDOODLE_LIKE` rendering defaults. Any newly introduced rendering override SHALL default to an explicit no-override state that follows the existing calculation exactly.

#### Scenario: Default visual rendering is unchanged
- **GIVEN** the default Current profile and an existing depiction fixture
- **WHEN** it is applied to production drawing items
- **THEN** the same drawing-style fields and existing fallback calculations are used as before the laboratory
- **AND** the production editor does not load Working profiles.

### Requirement: Working controls affect effective renderer properties
Every enabled Working control SHALL update a real `DrawingStyle` field or documented optional renderer override and SHALL refresh the actual Qt items immediately. A parameter that cannot be independently and safely calibrated SHALL be marked unavailable and SHALL NOT be represented by a decorative control.

#### Scenario: Edit a calibrated parameter
- **GIVEN** a visible fixture and a parameter supported by its renderer path
- **WHEN** its Working control changes
- **THEN** the Working item geometry/style changes on refresh, Current remains unchanged, and graph snapshots stay equal.

### Requirement: Profiles use a strict versioned schema
Profiles SHALL serialize typed dimensionless ratios relative to the reference bond length, discrete counts/options, schema and version. Loading SHALL reject invalid types, non-finite/out-of-range values and unknown fields. Compatibility SHALL be implemented through an explicit version migration registry. Profile files SHALL be stored/exported only in lab-owned locations and SHALL not overwrite application preferences.

#### Scenario: Invalid profile fails closed
- **GIVEN** JSON with unsupported version, unknown parameter, NaN/Infinity, wrong type, or out-of-range value
- **WHEN** the profile is loaded
- **THEN** the loader reports a structured validation error and does not partially apply it.

#### Scenario: Save and reopen Working profile
- **GIVEN** valid Working values
- **WHEN** the profile is saved and reopened
- **THEN** the parameter model is identical and renders deterministically for the same case/version.

### Requirement: Gallery fixtures are deterministic and chemically explicit
The gallery SHALL contain at least twelve stable case IDs covering labeled wedges, hashed/double junctions, multiple bonds/branches, rings and Haworth/aromatic depictions. Each fixture SHALL be identified as a chemically valid depiction case or a geometric-only fixture. Case construction SHALL not mutate graph coordinates after it is finalized.

#### Scenario: Select any gallery case
- **GIVEN** the gallery's complete case registry
- **WHEN** the user selects any stable case ID
- **THEN** the same deterministic graph and corresponding two-panel rendering are shown without manual setup.

### Requirement: View transforms and diagnostic overlays are presentation-only
Zoom and optional angle stress transforms SHALL be synchronized between panels and SHALL transform only the view. Diagnostic overlays SHALL be independent graphics items, removable, and excluded from ordinary exports.

#### Scenario: Rotate a scene for angular stress
- **GIVEN** a case with explicit stereo metadata
- **WHEN** the view is rotated through a supported angle
- **THEN** atom coordinates in both graphs and every stereo field remain unchanged.

### Requirement: Reproducible export records provenance
The lab SHALL export separate Current and Working PNGs, a before/after comparison, a gallery contact sheet, the exact JSON profile and a Markdown report containing stable case IDs, renderer/profile settings, ChemUSON version and source SHA. Export dimensions, backgrounds, framing and overlay inclusion SHALL be explicit.

#### Scenario: Repeat export
- **GIVEN** the same commit, profile, fixture, theme, zoom and export settings
- **WHEN** the same export command runs in the same controlled environment
- **THEN** images are deterministic within the documented Qt/font environment and report/profile data are byte-stable.

### Requirement: Lab lifecycle is bounded and isolated
The lab SHALL be launchable with a documented development command, and UI tests SHALL use offscreen Qt and explicit close/cleanup under external timeouts. Closing the lab SHALL not alter application state, documents, selections or undo history.

#### Scenario: Start and close offscreen
- **GIVEN** an offscreen Qt test process
- **WHEN** the window is constructed, exercised and closed
- **THEN** the process exits within its external timeout without a modal interaction, hang or abnormal signal.
