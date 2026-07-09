# clean-2d-baseline-report-dev-script Specification

## Purpose
TBD - created by archiving change add-clean2d-baseline-report-dev-script. Update Purpose after archive.
## Requirements
### Requirement: Developer Baseline Report Script
The system SHALL provide a developer-only Clean 2D baseline report script.

#### Scenario: Script is available
- **GIVEN** the project workspace
- **WHEN** a developer invokes the script
- **THEN** it provides commands for writing, comparing, and reviewing Clean 2D baseline reports.

### Requirement: JSON-Stable Report Generation
The script SHALL generate JSON-stable baseline reports.

#### Scenario: Write command creates report
- **GIVEN** an output path
- **WHEN** the write command runs
- **THEN** the script writes a JSON-stable baseline report containing schema, version, records, and summary.

### Requirement: Report Comparison
The script SHALL compare two baseline report JSON files using the existing report diff rules.

#### Scenario: Compare command runs
- **GIVEN** two baseline report JSON files
- **WHEN** the compare command runs
- **THEN** it prints stable diff JSON and indicates whether the reports are equivalent.

### Requirement: Review Policy Application
The script SHALL review two baseline report JSON files using the existing diff review policy.

#### Scenario: Review command runs
- **GIVEN** two baseline report JSON files
- **WHEN** the review command runs
- **THEN** it prints stable review summary JSON.

### Requirement: Explicit Exit Codes
The script SHALL use explicit exit codes for equivalent, changed and error outcomes.

#### Scenario: Command completes
- **GIVEN** the script command result
- **WHEN** the process exits
- **THEN** it returns `0` for success/equivalence, `1` for observable differences, and `2` for argument, read, write, or JSON errors.

### Requirement: Observational Outputs
The script SHALL keep all outputs observational-only in this change.

#### Scenario: Differences exist
- **GIVEN** compare or review detects differences
- **WHEN** output is printed
- **THEN** the output describes differences without enforcing geometry quality gates.

### Requirement: No Production Behaviour Changes
The script SHALL NOT change Clean 2D production behaviour, layout, ranking, selection, UI, backends or public document format.

#### Scenario: Script runs
- **GIVEN** the script executes report helpers
- **WHEN** Clean 2D is invoked through test-owned infrastructure
- **THEN** production behaviour is observed but not modified.

### Requirement: No Geometry Gates
The script SHALL NOT introduce geometry gates or CI enforcement in this change.

#### Scenario: Geometry-risk review item exists
- **GIVEN** review output contains geometry-risk
- **WHEN** the script exits with differences
- **THEN** the exit status only indicates observable differences, not failed geometry quality.

