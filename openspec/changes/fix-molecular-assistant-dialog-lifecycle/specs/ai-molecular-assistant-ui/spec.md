## ADDED Requirements

### Requirement: Dialog and per-job state SHALL follow Qt object lifetime

The Molecular Assistant UI SHALL clean its transient per-job state using a stable job ID, not by querying a dialog wrapper during teardown. For every job, cleanup SHALL be idempotent, abandon still-pending work, and remove the dialog, result, identity result, and transform context registries. A successful preview MAY retain its job state only while its dialog is alive and awaiting an explicit decision.

#### Scenario: Closing an unstarted dialog creates no stale job state
- **GIVEN** the Molecular Assistant dialog is open and no job has started
- **WHEN** the user closes it and Qt processes destruction
- **THEN** no synthetic job SHALL be created
- **AND** no dialog or per-job registry entry SHALL remain.

#### Scenario: Closing a pending dialog abandons by stable job ID
- **GIVEN** a job is associated with an open Molecular Assistant dialog
- **WHEN** the dialog is finished or destroyed while the worker is pending
- **THEN** the controller SHALL abandon that job by ID
- **AND** all four per-job registries SHALL remove that ID
- **AND** a later worker result SHALL NOT access the destroyed dialog or mutate a canvas.

#### Scenario: A completed preview remains live only until the dialog finishes
- **GIVEN** a worker has completed successfully and the dialog displays a preview
- **WHEN** the user closes/declines the dialog
- **THEN** result, identity, dialog, and transform context state SHALL be removed by ID
- **AND** no insertion SHALL occur.

#### Scenario: Accepted insertion cleans up without a stale wrapper
- **GIVEN** the user accepts a validated proposal and `WA_DeleteOnClose` destroys the dialog
- **WHEN** `finished` and `destroyed` callbacks run in either order
- **THEN** cleanup SHALL be idempotent and SHALL NOT query or call methods on a destroyed dialog.

#### Scenario: Repeated open and close leaves registries empty
- **GIVEN** the user opens and closes the assistant repeatedly, including dialogs with retries
- **WHEN** Qt processes all destruction events
- **THEN** all Molecular Assistant per-job registries SHALL be empty and no old job ID SHALL affect a new dialog.

#### Scenario: Transform dialog cleanup removes source context
- **GIVEN** a Molecular Assistant transform dialog has a job and source context
- **WHEN** that dialog closes or is destroyed
- **THEN** cleanup SHALL remove the transform context together with all other state for that job ID.

#### Scenario: Window shutdown tolerates a deleted dialog wrapper
- **GIVEN** ChemusonWindow begins shutdown with one or more Molecular Assistant dialogs and workers
- **WHEN** it abandons jobs, closes remaining live dialogs, and tracks worker termination
- **THEN** deleted dialog wrappers SHALL not be retained or called
- **AND** the window SHALL wait for owned workers under the existing shutdown protocol without a stale-dialog RuntimeError.

#### Scenario: A failure can be retried without retaining an old job registry
- **GIVEN** a provider failure is presented in a still-open dialog and the user retries
- **WHEN** a new worker job is started with a new ID
- **THEN** the prior job's registry state SHALL already be removed
- **AND** closing the dialog SHALL clean the new ID without consulting a mutable dialog job ID.
