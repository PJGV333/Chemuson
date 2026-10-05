## Why

The Molecular Assistant can propose and validate a new molecule, but it cannot transform an existing molecular structure. The user explicitly authorized Phase 5 as a whole-molecule transformation. This change implements the previously discussed minimal option: replace exactly one complete selected connected molecule, with review and explicit approval. It does not introduce per-atom editing operations or a non-destructive parallel product.

## What Changes

- Add a discoverable Structure-menu/command-palette action for transforming a selected molecule.
- Accept only a selection equal to one complete connected molecular component; reject partial or multiple-component selections.
- In the existing bounded worker, export the source graph to isolated ChemIO SMILES and pass that source plus the user's instruction through the existing M23 request contract.
- Review the source SMILES and validated proposed SMILES, then replace only after explicit approval.
- Keep source/document state unchanged on request failure, cancellation, decline, or stale source state; record replacement as one undoable operation.

## Capabilities

### New Capabilities
- `ai-whole-molecule-transformation`: Safe, reviewed transformation and undoable replacement of one complete selected molecule through existing M23 validation.

### Modified Capabilities
- `ai-molecular-assistant-ui`: Add the explicit whole-molecule transformation route while preserving the existing generate-and-insert flow.

## Impact

GUI action/command registry, Molecular Assistant dialog/controller adapter, main-window request context, canvas replacement orchestration, and offline tests. M23 request/result contracts and ChemIO validation remain unchanged. Clean2D remains diagnostic-only and is not called. No external dependency, persistence change, real provider call, agent/tool-calling behavior, or model-server requirement is introduced.
