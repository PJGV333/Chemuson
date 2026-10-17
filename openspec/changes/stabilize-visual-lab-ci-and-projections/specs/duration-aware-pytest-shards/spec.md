# Duration-aware pytest shards

## ADDED Requirements

### Requirement: Deterministic shards account for measured file duration
The CI planner SHALL partition the exact collected node-ID manifest by indivisible test file using a versioned per-file duration profile and a documented fallback for files without history. Given identical manifest, profile and shard count, it SHALL produce byte-equivalent assignment order and balanced estimated wall-time loads. It SHALL retain every collected node ID exactly once and MUST NOT skip/deselect tests.

#### Scenario: Balance known slow Qt and Clean2D files
- **GIVEN** test files with disparate measured durations per node ID
- **WHEN** the planner creates eight shards
- **THEN** high-cost files are distributed by estimated time rather than only node count
- **AND** assignments remain deterministic and coverage-exact.

### Requirement: Runtime evidence and bounded execution remain fail-closed
The shard runner SHALL record per-test/per-file elapsed time as tests complete and preserve partial timing evidence on timeout. It SHALL retain a finite external timeout, process-group termination, JUnit-to-assignment verification, and summary failure on any incomplete shard. Timing instrumentation SHALL NOT alter test outcomes or swallow process failures.

#### Scenario: Timeout leaves diagnostic evidence
- **GIVEN** a shard exceeds its configured timeout during a test run
- **WHEN** the runner terminates the process group
- **THEN** it records timeout, signal, last test and completed per-file durations
- **AND** CI summary remains failed because execution is incomplete.
