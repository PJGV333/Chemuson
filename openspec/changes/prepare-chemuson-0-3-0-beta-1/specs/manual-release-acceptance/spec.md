# Spec Delta — manual release acceptance

## Purpose

Provide a repeatable, versioned manual acceptance protocol for beta artifacts and stable promotion, with reproducible cases, explicit result states, severity and sign-off evidence.

## ADDED Requirements

### Requirement: Acceptance cases are reproducible and traceable
The release acceptance matrix SHALL assign stable case IDs and provide prerequisites/test data, step-by-step actions, expected results, a result field (`PASS`, `FAIL`, `BLOCKED`, or `NOT TESTED`), artifact version/SHA, tester, environment and evidence link. The same matrix SHALL be usable for beta and release-candidate artifacts.

#### Scenario: Tester records a case
- **WHEN** a tester executes a matrix case against an installed artifact
- **THEN** the recorded result identifies the exact artifact/version, environment and observed evidence
- **AND** unexecuted cases remain `NOT TESTED`, not implied PASS.

### Requirement: Linux AppImage acceptance distinguishes format from graphical behavior
The Linux distribution cases SHALL verify the AppImage Type 2 signature, successful `--appimage-extract`, expected AppDir entries and internal version before installation. Runner-side offscreen startup evidence SHALL be recorded separately from interactive visual acceptance; headless execution MUST NOT be reported as graphical QA.

#### Scenario: AppImage package is received
- **WHEN** a tester downloads a Linux preview/release package
- **THEN** its SHA/provenance and Type 2 structure are verified before installation
- **AND** visual, input, close/teardown and chemical workflows remain manual test cases.

### Requirement: Stable promotion has explicit severity gates
A stable promotion SHALL require approval of the owner and completion of the declared priority manual cases. Any open P0 or P1 attributable to the candidate SHALL block stable publication. P2/P3 issues MAY be accepted only with a documented workaround, owner decision and follow-up. The matrix SHALL distinguish product defects from known test-harness baseline exceptions.

#### Scenario: Critical defect is found
- **WHEN** a manual run reproduces data corruption, severe security exposure, launch failure, essential workflow failure, or serious chemical/connectivity/stereo mutation
- **THEN** the case is marked FAIL with P0/P1 severity
- **AND** stable promotion is blocked.

#### Scenario: Beta has an untested noncritical area
- **WHEN** an area is not yet tested on the beta artifact
- **THEN** it remains visibly `NOT TESTED` or `BLOCKED`
- **AND** beta availability may proceed only under the beta criteria, not by treating the missing evidence as a stable pass.

### Requirement: Chemistry and persistence acceptance checks preserve molecular data
Manual checks SHALL compare molecular connectivity, atom/bond properties, charges, aromaticity, stereo and coordinates across save/open, import/export, Undo/Redo and Clean2D as relevant. Clean2D acceptance SHALL reject silent chemical changes and SHALL not imply that all Clean2D work is complete.

#### Scenario: Persist and reopen a molecule
- **WHEN** a representative `.cmsn` document is saved, closed and reopened
- **THEN** its graph, chemistry-relevant properties and expected coordinates are preserved
- **AND** any mismatch is recorded as a failure rather than repaired silently.

### Requirement: Packaged icon defects block beta acceptance until manual retest
A confirmed P1 missing-icon defect SHALL be recorded as `FAILED — P1 blocks beta acceptance` per affected platform in the manual matrix. For this campaign it MUST block beta publication, not only stable promotion, until the owner completes the manual retest. Automated source/offscreen or frozen-resource smoke checks MUST NOT turn those manual cases into PASS. The owner MUST install the corrected Windows portable and Linux portable/AppImage artifacts and visually verify pointer/select, single bond, aromatic ring, search, undo, redo, new document and clean tools in light and dark themes before changing their status.

#### Scenario: Frozen smoke passes while manual retest is pending
- **WHEN** the packaged executable reports all required SVGs and visible QtSvg rasters
- **THEN** this is recorded only as automated packaging evidence
- **AND** the affected manual cases remain blocked/failed until owner visual retest on the real packages.

#### Scenario: Owner confirms icons on corrected packages
- **WHEN** the owner records artifact SHA, platform, themes and screenshots/evidence after retest
- **THEN** only the corresponding platform/theme matrix case may change from FAILED to PASS.

### Requirement: Packaged ChemName regression stays blocked until real package retest
The manual acceptance matrix SHALL record the frozen-package ChemName template omission as a P1, distinguish Windows/Linux PyInstaller evidence from the Flatpak package, and keep corrected package retests `NOT TESTED` until the owner evaluates the exact preview SHAs. Automated source/frozen name checks MUST NOT be reported as GUI/manual acceptance.

#### Scenario: Corrected ChemName package passes automation
- **WHEN** the new Windows/Linux executable and Flatpak smoke checks pass
- **THEN** their names/resources are recorded as automated evidence only
- **AND** the owner's manual status-bar, annotation, and package retest remains pending.

### Requirement: AI acceptance distinguishes structure validity from identity
Molecular Assistant cases SHALL record provider availability, source provenance, ChemIO validity and identity/reference status separately. A valid SMILES alone SHALL NOT be considered proof of the requested molecular identity. Tests SHALL cover offline/provider failure behavior without assuming a Qwen endpoint exists.

#### Scenario: AI proposal and reference differ
- **WHEN** the assistant proposes a valid structure with a different molecular identity from the reference
- **THEN** the tester records both provenance and mismatch behavior
- **AND** the candidate is not reported as identity-verified solely because ChemIO accepted it.

#### Scenario: AI provider is unavailable
- **WHEN** a model endpoint is offline or returns invalid/exhausted output
- **THEN** the application fails controllably or offers an explicitly identified valid reference
- **AND** the test records the observed result without claiming a model success.
