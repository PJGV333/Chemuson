# clean-2d-campaign1-benchmark-observability Specification

## Purpose

Campaign 1 defines the reproducible, test-owned benchmark and observability foundation required before algorithmic Clean2D changes.

## ADDED Requirements

### Requirement: Stable benchmark case identity and taxonomy

The benchmark SHALL preserve each regression case's stable name and SHALL expose an explicit size class, family, and classification tags without deriving execution behavior from a case name.

#### Scenario: Case metadata is stable

- **GIVEN** the regression corpus is enumerated twice
- **WHEN** case metadata is read
- **THEN** names, size classes, families, and tags are identical and names remain unique

#### Scenario: Existing case IDs remain visible

- **GIVEN** an existing known-delicate, known-failure, stereo-sensitive, selection-boundary, or complex-policy case
- **WHEN** Campaign 1 metadata is emitted
- **THEN** its existing name and classification tags remain present

### Requirement: Derived topology metadata is observational

The benchmark SHALL derive measurable topology metadata from the built graph and SHALL represent unavailable values explicitly rather than guessing them.

#### Scenario: Graph metadata is reproducible

- **GIVEN** a valid corpus case builder
- **WHEN** its graph metadata is derived twice
- **THEN** atom count, heavy-atom count, bond count, ring count, and connected-component count are equal

#### Scenario: Metadata does not route Clean2D

- **GIVEN** a corpus case with a family or topology label
- **WHEN** Clean2D executes
- **THEN** the metadata SHALL be observational and SHALL NOT select a candidate or alter layout behavior

### Requirement: Baseline records expose one diagnostic metric vector

Each baseline record SHALL expose the existing geometry metrics together with candidate count, candidate sources, result state, stable reason, runtime evidence, and topology metadata using JSON-safe values.

#### Scenario: Complete evidence record

- **GIVEN** a corpus case executes through the existing Clean2D engine
- **WHEN** a baseline record is built
- **THEN** it contains the case identity, taxonomy, topology metadata, before/after diagnostic metrics, candidate evidence, result state, stable reason, and runtime evidence

#### Scenario: Unavailable metric is explicit

- **GIVEN** a metric cannot be computed for a case
- **WHEN** the evidence record is serialized
- **THEN** the metric is represented as `null` and no NaN or infinity is emitted

### Requirement: Baseline execution is deterministic apart from runtime

Repeated baseline execution SHALL produce equivalent canonical records apart from explicitly ephemeral runtime evidence.

#### Scenario: Repeated reports compare equivalent

- **GIVEN** the same source tree, corpus, mode, parameters, and seed
- **WHEN** the baseline runner writes two reports
- **THEN** comparison SHALL report equivalent unless a non-runtime observable field changed

#### Scenario: Runtime remains evidence only

- **GIVEN** two runs have different wall-clock durations
- **WHEN** their reports are compared
- **THEN** the runtime difference SHALL be retained as evidence but SHALL NOT alone make reports non-equivalent

### Requirement: Baseline CLI has auditable outcomes

The developer baseline CLI SHALL provide stable machine-readable write, compare, and review operations with distinct success, changed, and error exit codes.

#### Scenario: Write report

- **GIVEN** a valid output path
- **WHEN** the developer invokes `write`
- **THEN** a versioned JSON report is written with one record per corpus case

#### Scenario: Compare reports

- **GIVEN** two valid reports
- **WHEN** the developer invokes `compare`
- **THEN** JSON diff output identifies added, removed, or changed cases and the exit code distinguishes equivalent from changed reports

#### Scenario: Review changed reports

- **GIVEN** two reports with observable differences
- **WHEN** the developer invokes `review`
- **THEN** the output classifies the changed fields and does not silently promote the change

### Requirement: Campaign 1 does not alter production behavior

Campaign 1 SHALL remain test-owned and observational.

#### Scenario: Production behavior remains outside scope

- **GIVEN** Campaign 1 benchmark helpers and reports are not imported by normal application execution
- **WHEN** Clean2D runs outside the test runner
- **THEN** candidate generation, ranking, layout, GUI behavior, persistence, and chemical data remain unchanged
