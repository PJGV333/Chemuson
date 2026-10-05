# GUI Asynchronous Worker Shutdown Specification

## Purpose

Specify safe, deferred shutdown of asynchronous Qt workers owned by `ChemusonWindow` or its child canvases/controllers without changing worker computation semantics.

## ADDED Requirements

### Requirement: Close approval precedes worker shutdown

`ChemusonWindow` SHALL preserve its existing dirty-document and pending-update close decisions before entering worker shutdown. If any decision rejects close, active jobs SHALL remain usable and SHALL NOT be abandoned merely because close was attempted.

#### Scenario: User cancels dirty-document close
- **WHEN** at least one document is dirty and the user cancels its close confirmation
- **THEN** the close event is ignored
- **AND** asynchronous jobs continue normally
- **AND** no shutdown flag, interruption request, or job abandonment is applied.

#### Scenario: All normal close gates approve
- **WHEN** every dirty-document and pending-update close gate approves
- **THEN** the window enters shutdown exactly once
- **AND** no new window-owned background jobs start.

### Requirement: Active QThreads outlive their QObject owners

The window SHALL NOT be destroyed while any QThread in its QObject descendant tree is running. The window MAY defer close until all such threads finish naturally. It MUST NOT use `QThread.terminate()` or a busy loop. `requestInterruption()` SHALL be treated only as a cooperative request and SHALL NOT be represented as cancelling blocking work that does not check interruption.

#### Scenario: A descriptor worker is active
- **WHEN** close is approved while `DescriptorWorker` is executing isolated ChemIO/RDKit work
- **THEN** the pending properties timer cannot start another descriptor job
- **AND** the window remains alive until the QThread stops
- **AND** the thread/controller/window are destroyed only after thread completion.

#### Scenario: A bounded blocking worker is active
- **WHEN** Name→Structure, Molecular Assistant, CompChem, template SMILES export, or canvas analysis is active during approved close
- **THEN** shutdown may request cooperative interruption
- **AND** blocking work may finish under its existing finite timeout
- **AND** close is deferred until its QThread has stopped.

### Requirement: Late worker results are inert during shutdown

Once shutdown begins, worker results SHALL NOT update the canvas, document, selection, undo stack, docks, status bar, dialogs, or progress UI. Shutdown SHALL close/discard progress affordances safely and SHALL suppress late errors/confirmation dialogs as well as successful results.

#### Scenario: Name→Structure completes late
- **WHEN** a Name→Structure worker returns after shutdown begins
- **THEN** no QMessageBox is shown
- **AND** no graph is inserted or canvas/status widget updated
- **AND** its thread references are safely retired after completion.

#### Scenario: Molecular Assistant or CompChem completes late
- **WHEN** an assistant or CompChem worker emits a result/frame after shutdown begins
- **THEN** the result/frame is not relayed into the UI or document
- **AND** the QThread still completes and is retained until safe destruction.

### Requirement: Deferred close completes without repeating user decisions

If close is approved while workers are active, the original close event SHALL be ignored until completion. After all worker QThreads stop, the window SHALL finish closing without prompting again, restarting jobs, or leaving invalid Qt references.

#### Scenario: Final worker completes
- **WHEN** the last active descendant QThread has finished
- **THEN** the window completes the approved close
- **AND** no warning `QThread: Destroyed while thread is still running` or process crash occurs.

#### Scenario: No workers are active
- **WHEN** close is approved and no descendant QThread is running
- **THEN** the window closes through the existing path without unnecessary deferral.

### Requirement: Lifecycle regression tests are deterministic and offline

Tests SHALL use fake/slow finite workers and SHALL NOT require live providers, network access, or model servers. Coverage SHALL include each window-owned worker family, cancelled close, empty-worker close, and the original order-dependent crash pair in both orders.

#### Scenario: Original problematic test order
- **WHEN** the CompChem export/window teardown test is followed by the Molecular Assistant worker test
- **THEN** both pass repeatedly and the process exits cleanly.

#### Scenario: Reverse order
- **WHEN** the two tests run in reverse order
- **THEN** both continue to pass.
