## ADDED Requirements

### Requirement: smart-repair-before-propose
The system SHALL check if the base geometry is distorted before proposing alternatives in the `smart propose` flow. If a safe repair exists, the system SHALL apply it first.

#### Scenario: Repair applied to distorted geometry
- **WHEN** the user initiates a `smart propose` on a distorted geometry with a known safe repair
- **THEN** the system SHALL apply the repair to the base geometry before generating any alternatives

## MODIFIED Requirements

### Requirement: clean-2d-smart-propose
The `smart propose` flow SHALL ensure that the base geometry is in a valid state (repaired if necessary) before proposing alternatives to ensure bond length requirements are met.

#### Scenario: Valid bond lengths in proposals
- **WHEN** `smart propose` is executed on a distorted geometry
- **THEN** the proposed alternatives SHALL have bonds that meet the minimum length requirements because the base geometry was repaired first
