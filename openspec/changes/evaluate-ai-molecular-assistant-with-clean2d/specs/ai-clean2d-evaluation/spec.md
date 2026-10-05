# AI Clean2D Evaluation Specification

## Purpose

Define an explicit, non-mutating campaign evaluator that passes only ChemIO-validated M23 proposals through the existing Clean2D engine and reports before/after visual metrics without changing Clean2D behavior.

## ADDED Requirements

### Requirement: Explicit Network Invocation
The evaluator SHALL make provider requests only when its CLI is explicitly invoked with a non-empty description, endpoint, and model. Importing the evaluator or using the normal ChemUSON application SHALL NOT initiate provider I/O.

#### Scenario: Tool is imported
- **GIVEN** the evaluator module is imported
- **WHEN** no command is invoked
- **THEN** no provider, transport, network request, or model service is started.

#### Scenario: Configuration is incomplete
- **GIVEN** endpoint, model, or description is missing or invalid
- **WHEN** the CLI arguments are validated
- **THEN** evaluation is not started and no HTTP request is made.

### Requirement: Validated-Graph Gate
The evaluator SHALL invoke Clean2D only for an M23 result whose status is `success`, whose `validation_passed` value is true, and whose graph is a non-empty `MolGraph`.

#### Scenario: M23 returns a failure
- **GIVEN** generation, decoding, or isolated ChemIO validation fails
- **WHEN** the evaluator builds its report
- **THEN** the report preserves the stable M23 status/reason, contains no Clean2D metrics, and the Clean2D engine is not called.

#### Scenario: M23 returns a validated graph
- **GIVEN** M23 returns a successful, validated non-empty `MolGraph`
- **WHEN** evaluation proceeds
- **THEN** the graph is passed to the existing M02 Clean2D engine without invoking Clean2D from M23.

### Requirement: Before-and-After Diagnostic Metrics
The evaluator SHALL measure the existing Clean2D layout-quality metrics on the validated graph's initial coordinates and, when an accepted result exists, on the selected candidate coordinates. The report SHALL preserve Clean2D's returned result state and stable reason independently of metric values.

#### Scenario: Accepted candidate
- **GIVEN** Clean2D returns an accepted selected candidate
- **WHEN** the evaluation report is produced
- **THEN** it includes `before` and `after` records with the established quality class, reason, crossings, minimum non-bonded distance, minimum ring degeneracy, normalized bond-length RMS/max error, angle RMS/max deviation, and visual score.

#### Scenario: No accepted candidate
- **GIVEN** Clean2D returns no accepted candidate
- **WHEN** the report is produced
- **THEN** `before` metrics are present, `after` is null, and the engine's state/reason are preserved.

#### Scenario: Metrics appear worse
- **GIVEN** diagnostic metrics are poor or worsen after evaluation
- **WHEN** the report is serialized
- **THEN** the metrics do not alter the Clean2D state, candidate selection, or report success semantics.

### Requirement: Diagnostic-Only Non-Mutation
The evaluator SHALL NOT mutate the validated `MolGraph` or editor/document state. It SHALL NOT insert a graph into the canvas.

#### Scenario: Clean2D evaluation completes
- **GIVEN** a validated graph is evaluated
- **WHEN** candidate coordinates and metrics are collected
- **THEN** source graph chemistry, coordinates, and M23 result remain unchanged and no GUI/canvas API is called.

### Requirement: JSON-Safe Private Report
The evaluator SHALL emit one JSON-safe report to stdout and SHALL NOT include the natural-language description, API key, full HTTP diagnostics, or unparsed provider response.

#### Scenario: Report contains unavailable geometry values
- **GIVEN** a metric is non-finite or not applicable
- **WHEN** the report is encoded as JSON with NaN disallowed
- **THEN** the metric is null and serialization succeeds.

#### Scenario: Provider fails
- **GIVEN** M23 returns a controlled provider failure
- **WHEN** the report is serialized
- **THEN** only stable status/reason and non-secret provenance fields are exposed.

### Requirement: Offline Verification
Tests SHALL use fake providers/transports and injected engine behavior. Tests SHALL NOT require external APIs, credentials, network access, or a manually started model service.

#### Scenario: Test suite executes
- **GIVEN** focused evaluator tests run in CI or locally
- **WHEN** provider and engine dependencies are exercised
- **THEN** deterministic fakes cover success, failure, call ordering, non-mutation, and finite JSON serialization.
