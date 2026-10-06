## Why

The Molecular Assistant can propose and validate a new molecule, but it cannot transform an existing molecular structure. The user explicitly authorized Phase 5 as a whole-molecule transformation. This change implements the previously discussed whole-molecule transformation, with review and explicit approval. Preview offers either inserting a separate variant while preserving the source or replacing exactly one complete selected connected molecule. It does not introduce per-atom editing operations.

## What Changes

- Add a discoverable Structure-menu/command-palette action for transforming a selected molecule.
- Accept only a selection equal to one complete connected molecular component; reject partial or multiple-component selections.
- In the existing bounded worker, export the source graph to isolated ChemIO SMILES and pass a typed `MolecularTransformationRequest(source_smiles, instruction)` through M23's provider-neutral transform operation.
- Review the source SMILES and validated proposed SMILES, then explicitly insert a separate variant or replace the source.
- Keep source/document state unchanged on request failure, cancellation, decline, or stale source state. Variant insertion preserves the source; replacement and variant insertion are each one undoable operation.

## Capabilities

### New Capabilities
- `ai-whole-molecule-transformation`: Safe, reviewed transformation and undoable replacement of one complete selected molecule through existing M23 validation.

### Modified Capabilities
- `ai-molecular-assistant-ui`: Add the explicit whole-molecule transformation route while preserving the existing generate-and-insert flow.

## Impact

GUI action/command registry, Molecular Assistant dialog/controller adapter, main-window request context, canvas replacement orchestration, and offline tests. M23 gains a typed provider-neutral transformation request while reusing the existing strict decoder and ChemIO validation; no duplicate validation path is added. Clean2D remains diagnostic-only and is not called. No external dependency, persistence change, real provider call, agent/tool-calling behavior, or model-server requirement is introduced.
