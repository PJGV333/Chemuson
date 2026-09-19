## ADDED Requirements

### Requirement: Operational Resilience Ownership Is Explicit

The repository SHALL document M22 as the owner of generic crash logging,
autosave, recovery snapshots and failure isolation; M08/M10 as owners of Qt
task containment; M14 as the owner of update telemetry; and M19 as the bootstrap
installer of resilience hooks.

#### Scenario: Ownership audit

- **WHEN** the operational resilience audit is inspected
- **THEN** each concern has exactly one documented primary owner

### Requirement: No Premature Module Is Created

The audit SHALL NOT create a new placeholder module for cross-cutting
resilience concerns while the existing owners remain cohesive.

#### Scenario: Reserved slots

- **WHEN** the catalog is inspected after the audit
- **THEN** M23 and M24 remain reserved and no new operational-resilience package
  exists
