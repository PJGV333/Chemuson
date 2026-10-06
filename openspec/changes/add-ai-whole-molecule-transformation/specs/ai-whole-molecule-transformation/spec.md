# AI Whole-Molecule Transformation Specification

## Purpose

Specify a reviewed, undoable transformation of exactly one complete selected connected molecule using the existing M23 generation and ChemIO validation contract.

## ADDED Requirements

### Requirement: A transform request SHALL target one complete selected component

The GUI SHALL start transformation only when the selected atom IDs are exactly one complete connected component in the active canvas. It SHALL NOT infer or silently expand a partial selection. Empty selection, multiple selected components, or selected bonds outside the component SHALL be rejected without provider I/O.

#### Scenario: Complete connected component selected
- **WHEN** the user invokes transformation with exactly one complete connected component selected
- **THEN** the request is enabled for that source component
- **AND** the source graph is snapshotted before asynchronous work begins.

#### Scenario: Selection is empty, partial, or ambiguous
- **WHEN** no atoms are selected, only part of a connected molecule is selected, multiple components are selected, or selection IDs are stale
- **THEN** the GUI reports that one complete molecule must be selected
- **AND** no worker/provider request starts
- **AND** the document remains unchanged.

### Requirement: Source export and M23 generation SHALL be asynchronous and bounded

The GUI SHALL export the selected source graph to isolated ChemIO SMILES and compose that complete source SMILES with the user's transformation instruction inside the existing Molecular Assistant worker. The existing provider timeout, M23 strict response decoding, isolated validation, result types, and size limits SHALL remain in force. The user-visible source/proposed molecules and API key SHALL NOT be persisted.

#### Scenario: Source SMILES export succeeds
- **WHEN** the worker exports the source molecule
- **THEN** the same worker sends the composed instruction through the existing M23 generation/validation path
- **AND** the GUI thread remains responsive.

#### Scenario: Source SMILES export fails
- **WHEN** isolated source export fails or times out
- **THEN** provider generation is not called
- **AND** a controlled failure is shown without raw subprocess details
- **AND** no canvas/document mutation occurs.

### Requirement: The proposal SHALL be reviewed before replacement

The review UI SHALL show source SMILES, proposed SMILES, provider/model provenance when available, and the existing caveat that ChemIO parser acceptance does not prove semantic correctness. It SHALL offer Insert Variant, Replace Original, and Cancel. Mutation SHALL occur only after an explicit user choice.

#### Scenario: User declines or closes the preview
- **WHEN** a user declines the transformation or closes the dialog
- **THEN** the source molecule, other graph components, selection, undo index, and dirty state remain unchanged
- **AND** any late result is ignored.

#### Scenario: User chooses Replace Original
- **WHEN** the user explicitly chooses Replace for a successful, ChemIO-validated proposal
- **THEN** only the selected source component is removed and the proposal is placed at the source component's original center
- **AND** the replacement is recorded as exactly one undoable operation.

#### Scenario: User chooses Insert Variant
- **WHEN** the user explicitly chooses Insert Variant for a successful, ChemIO-validated proposal
- **THEN** the complete source component remains unchanged
- **AND** the proposal is inserted as a separate molecule at a non-overlapping position relative to the source
- **AND** the variant insertion is exactly one undoable operation without automatic Clean2D.

### Requirement: Replacement SHALL be conditional on an unchanged source

Before applying an approved proposal, the application SHALL verify that the target canvas remains open and active, the same complete component is selected, and all source Atom/Bond values and coordinates still match the request snapshot. If any check fails, replacement SHALL be rejected without mutation.

#### Scenario: Source document or selection changed while generation was pending
- **WHEN** the user changes the source molecule, selection, or active document before approval
- **THEN** the application shows a controlled stale-source notice
- **AND** no atoms or bonds are deleted or inserted.

### Requirement: Undo and redo SHALL restore the whole-molecule transaction

The replacement SHALL compose existing delete/add canvas commands in one undo macro. Undo SHALL restore the exact source component and remove the proposal after Replace; for Insert Variant, Undo SHALL remove only the variant. Redo SHALL replay the selected operation. Other molecules and unrelated drawing objects SHALL remain untouched.

#### Scenario: User undoes and redoes replacement
- **WHEN** a successful replacement is followed by Undo and Redo
- **THEN** Undo/Redo restores the original source for Replace, or removes/restores only the inserted variant for Insert Variant
- **AND** the source is untouched by variant Undo/Redo
- **AND** each chosen operation occupies one undo-stack step.

### Requirement: Tests SHALL remain deterministic and offline

Tests SHALL inject fake SMILES export and generation behavior. They SHALL NOT invoke real providers, network services, model servers, or manual model sessions.

#### Scenario: Offline transform tests execute
- **WHEN** focused transformation tests run
- **THEN** request composition, failures, stale-state rejection, explicit approval, and undo/redo are deterministic without external I/O.
