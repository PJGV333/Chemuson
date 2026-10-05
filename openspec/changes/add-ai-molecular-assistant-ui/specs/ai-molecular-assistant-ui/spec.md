## ADDED Requirements

### Requirement: Structure generation SHALL be discoverable from the existing UI

ChemUSON SHALL expose a minimal AI structure-generation action in the existing Structure menu and command palette. The action SHALL use the existing GUI composition and SHALL NOT add a chat interface.

#### Scenario: User opens AI structure generation
- **WHEN** the user selects the AI structure action from the Structure menu or command palette
- **THEN** ChemUSON opens the same minimal request flow
- **AND** the existing drawing, document, and Clean2D actions remain available

### Requirement: Requests SHALL use explicit, transient provider configuration

The UI SHALL collect a non-empty molecular description and the explicit endpoint/model configuration required by M23's OpenAI-compatible adapter. An optional API key SHALL be masked while entered and SHALL remain transient to the active request; the UI SHALL NOT persist prompts, API keys, or provider responses.

#### Scenario: Required request configuration is missing
- **WHEN** the user submits an empty description, endpoint, or model
- **THEN** generation SHALL NOT start
- **AND** the UI SHALL identify the missing required value without contacting a provider

#### Scenario: API key is entered
- **WHEN** the user enters an optional API key
- **THEN** the input is masked
- **AND** the key is passed only to the current M23 provider configuration and is not written to settings, logs, or the document

### Requirement: Provider generation and ChemIO validation SHALL not block the GUI

The UI SHALL run M23 generation and its isolated ChemIO validation away from the GUI thread, using the configured finite provider timeout. The UI SHALL remain responsive while a request is pending, and a late result for a closed or cancelled UI request SHALL NOT be applied to a canvas.

#### Scenario: Generation is pending
- **WHEN** a request is sent to the provider
- **THEN** the request and validation execute outside the GUI thread
- **AND** the GUI event loop remains responsive
- **AND** the user can close the request flow without inserting a late result

#### Scenario: Provider or validation fails
- **WHEN** M23 returns any non-success result or the worker reports an unexpected failure
- **THEN** the UI displays a controlled status/reason without exposing raw transport or exception diagnostics
- **AND** it SHALL NOT call molecule insertion

### Requirement: A validated proposal SHALL be reviewed before insertion

The UI SHALL present a successful result for review before any document mutation. The review SHALL identify provider/model when available, show the exact proposed SMILES, indicate ChemIO parser acceptance, and state that parser acceptance does not prove semantic correctness. The UI SHALL insert only after an explicit user choice.

#### Scenario: User declines a successful proposal
- **WHEN** the user closes or declines the result preview
- **THEN** no graph is inserted
- **AND** the active canvas, selection, undo history, and dirty state remain unchanged

#### Scenario: User approves a successful proposal
- **WHEN** the user explicitly chooses Insert on a successful preview
- **THEN** the validated `MolGraph` is inserted into the active canvas through the existing normal insertion path
- **AND** the inserted atoms are selected
- **AND** the insertion is recorded as one undoable operation

### Requirement: Failed or abandoned generation SHALL preserve editor state

Any invalid request, provider failure, malformed response, invalid structure, validation failure, cancellation, late result, or user decline SHALL leave the active canvas graph, selection, undo index, and dirty/clean state unchanged.

#### Scenario: Generation fails before preview
- **WHEN** a request does not produce a successful validated graph
- **THEN** the canvas graph and selection are unchanged
- **AND** the undo index and dirty/clean state are unchanged

#### Scenario: Undo and redo operate on an approved insertion
- **WHEN** the user inserts a reviewed proposal and then invokes Undo and Redo
- **THEN** Undo removes the inserted molecular graph and restores the prior clean/dirty state
- **AND** Redo restores the inserted graph as one normal editor operation
- **AND** selection follows the canvas's existing undo/redo behavior without retaining references to removed items

### Requirement: UI integration tests SHALL remain offline and deterministic

The GUI integration SHALL be testable using fake providers/services and SHALL NOT require network access, model servers, credentials, or manual interaction.

#### Scenario: A fake successful assistant result is reviewed before dispatch
- **WHEN** an offline test supplies a successful M23 result
- **THEN** the preview is observable without mutation
- **AND** explicit Insert dispatches the validated graph to the normal canvas insertion API with inserted-item selection enabled
- **AND** the canvas's existing undo/redo behavior remains the insertion mechanism

#### Scenario: A fake failure result is handled
- **WHEN** an offline test supplies each relevant failure class or abandons a pending job
- **THEN** the controlled UI outcome is observable and the editor-state snapshot is unchanged