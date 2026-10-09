# Delta de especificación — pytest CI reproducible

## Purpose

Definir contratos verificables para las regresiones de pytest conocidas y la ejecución completa de CI mediante grupos acotados con cobertura explícita.

## ADDED Requirements

### Requirement: Clean2D regressions verify backend invocation independently of candidate deduplication

The cyclic Clean2D test SHALL prove that the isolated RDKit backend is invoked and that its candidate geometry is valid even when geometric deduplication retains an equivalent candidate from another source. The test MUST NOT infer backend non-invocation from the absence of a source label after deduplication.

#### Scenario: Cyclic template duplicates isolated RDKit layout
- **GIVEN** a cyclic graph whose isolated RDKit layout is geometrically equivalent to an internal template
- **WHEN** publication candidates are generated
- **THEN** the real isolated backend is observed exactly once and yields a valid candidate
- **AND** the returned collection may retain only one candidate for the shared geometry.

### Requirement: Qt preferences tests isolate the actual QSettings route

Tests whose expected preferences are defaults SHALL isolate the NativeFormat/UserScope path after QApplication creation and restore the previous path regardless of test outcome. Changing `XDG_CONFIG_HOME` alone SHALL NOT be treated as sufficient when Qt may have cached standard paths.

#### Scenario: Persisted preference exists outside the test
- **GIVEN** a user preference such as `resolution_method=ai` exists in the previously selected NativeFormat path
- **WHEN** the isolated dialog test starts
- **THEN** it observes its fresh documented defaults
- **AND** the external preference and route are restored unchanged after the test.

### Requirement: Packaged-worker smoke probes are collection-order independent

Checks that assert the absence/presence of Qt application hooks or `chemuson.gui` modules SHALL execute in fresh subprocesses. The real package smoke SHALL continue to fail closed when a GUI module or active GUI application is present.

#### Scenario: Another test imports GUI modules during collection
- **GIVEN** pytest has imported GUI modules in its parent interpreter
- **WHEN** the packaged smoke import contracts run
- **THEN** each probe evaluates only its own fresh interpreter state
- **AND** active QApplication or `chemuson.gui` remains detectable within that interpreter.

### Requirement: CompChem async lifecycle is observable without increasing its timeout

The async CompChem regression SHALL observe `worker.finished`, `QThread.finished`, and controller `job_finished` with the existing bounded wait. On failure, the observed event order and active job state SHALL be available in the assertion report.

#### Scenario: Worker result and thread finish race
- **WHEN** a fake backend completes asynchronously
- **THEN** the result, worker signal, thread completion, and controller completion are observed in a coherent order within the existing bound
- **AND** a missing signal is reported as a failure, not converted into success by a longer timeout.

### Requirement: Side-panel layout assertions measure stable geometry

Primary-tab assertions SHALL be evaluated only after consecutive event-loop observations show stable geometry. The test SHALL retain the existing minimum font, padding, gap, width and unclipped-viewport requirements.

#### Scenario: Tab appears outside the viewport during layout
- **WHEN** window geometry changes
- **THEN** the test waits a finite bounded interval for stable geometry and records scroll/viewport measurements
- **AND** it passes only if every tab remains visible, padded, separated and unclipped in stable geometry.

### Requirement: Every collected pytest node executes exactly once in bounded CI shards

CI SHALL compare the complete pytest collection against a versioned sorted node-ID manifest, partition that exact set deterministically, execute every node exactly once across eight shards, and verify each shard's JUnit count against its assignment. Every shard SHALL have an external 300-second process limit and produce a machine-readable result. CI MUST fail on collection mismatch, duplicate/missing node IDs, timeout, signal/abort, missing artifacts, failed tests or invalid result counts. It MUST NOT use exclusions, skips introduced to hide failures, or `continue-on-error`.

#### Scenario: Complete shard campaign succeeds
- **GIVEN** the repository collects the same node IDs as the committed manifest
- **WHEN** CI plans and executes the shard matrix
- **THEN** each collected node appears in exactly one shard assignment and one JUnit execution
- **AND** all eight reports validate and every existing platform smoke job remains required.

#### Scenario: A pytest process aborts or exceeds its limit
- **WHEN** a shard times out, exits by signal, fails, or omits its report
- **THEN** the report/summary names the shard and available last-test/log information
- **AND** the workflow fails without suppressing the original outcome.

## Invariants

- No chemistry, GUI, preference, worker-production, release, Windows-smoke or Flatpak-smoke behavior changes.
- No test ID is excluded or omitted; each collected node is assigned once.
- External shard timeout is at most 300 seconds.
- A stable negative tab position remains a test failure until a narrowly justified visual correction is independently demonstrated.
