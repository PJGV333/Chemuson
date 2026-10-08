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
