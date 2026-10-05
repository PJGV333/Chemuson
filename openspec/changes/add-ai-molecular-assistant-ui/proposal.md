## Why

M23 can already turn a natural-language description into a ChemIO-validated `MolGraph`, but users cannot invoke it from ChemUSON. A small, explicit UI path is needed to review the proposal before it enters a document, without coupling the assistant to canvas mutation or Clean2D.

## What Changes

- Add a discoverable Structure action and command-palette entry for requesting a molecular structure.
- Collect the description and explicit OpenAI-compatible endpoint/model options for the current request; mask any API key and do not persist it.
- Run generation and isolated validation off the GUI thread, then present provider/model, SMILES, validation state, and the semantic-correctness caveat before insertion.
- Insert only after explicit user approval, through the canvas's existing undoable molecule insertion path.
- Add deterministic UI/controller tests for success, rejection, cancellation, undo/redo, and preservation of document state on failure.

## Capabilities

### New Capabilities
- `ai-molecular-assistant-ui`: Minimal GUI invocation, asynchronous result review, and explicit undoable insertion of an M23-validated molecular graph.

### Modified Capabilities
- None. M23's provider/validation contract and Clean2D requirements remain unchanged.

## Impact

Likely affected areas: M08 GUI actions, command registry, main-window orchestration, a small GUI controller/worker, and `architecture/modules.yml` to record M08's new M23 dependency. Tests will use fake providers/transports and must not access a live endpoint. No external dependency, persistent credential setting, chat history, Clean2D behavior, or M23 API change is introduced.