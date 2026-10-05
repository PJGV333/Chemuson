## Context

M23 accepts one natural-language description and returns a `MolGraph` only after its strict SMILES decoder and isolated ChemIO validation succeed. Its contract has no document or canvas context. The GUI already has a modeless worker/controller flow, isolated ChemIO SMILES export, `DeleteSelectionCommand`, and undoable canvas insertion helpers. Phase 4.5's Clean2D evaluation is a standalone diagnostic and does not authorize calling Clean2D during this transformation.

## Goals / Non-Goals

**Goals:**
- Transform exactly one selected whole connected component, never a partial fragment or multiple disconnected molecules.
- Keep source SMILES export, provider generation, and M23 validation off the GUI thread, under their existing finite timeouts.
- Show both source and validated proposal before the user explicitly chooses Replace.
- Reject replacement if the source component or selection changed while the request was pending.
- Replace the component as one undo/redo step while preserving other molecules and canvas objects.
- Test every path with fakes; no live provider or network service.

**Non-Goals:**
- Atom/bond operation DSLs, AI-authored code/tools, automatic replacement, chat/history, persistence of prompts/keys, or provider contract changes.
- Transformation of partial selections, multiple components, or an entire drawing containing separate molecules.
- Clean2D/geometry optimization or any change to chemical validation semantics.

## Decisions

1. **Interpret the user authorization as the existing option 1.** The UI accepts exactly one full connected selected component. The selected atom-ID set must equal that component. Partial selections, empty selection, multiple components, and out-of-component selected bonds are rejected before a worker starts. Do not guess omitted atoms or expand a partial selection silently.
2. **Reuse the existing M23 contract.** A per-job request transform runs in the Molecular Assistant worker: export a deep copy of the selected source graph with `molgraph_to_smiles_isolated_or_error(..., timeout_s=8.0)`, then compose a bounded instruction containing that complete source SMILES and the user's transformation instruction. The existing generator, provider timeout, strict response decoder, isolated ChemIO validation, result type, and finite limits remain authoritative. Source-export errors fail closed before the provider call and do not expose raw exception data.
3. **Keep the ordinary generation flow unchanged.** Add a separate discoverable QAction that reuses the same dialog/controller and provides a per-job description-transform hook. Ordinary draw-new requests pass no hook and preserve existing behavior. The API key remains transient and is never included in the transformation context or report/logs.
4. **Snapshot and revalidate before destructive replacement.** Capture the target canvas, exact selected IDs, complete source graph signature (all Atom/Bond dataclass fields, including coordinates), and source bounding-box center. Before Replace, require that the target document remains open and active, the same full component remains selected, and its complete signature still matches. A stale/missing target produces a controlled notice and no mutation.
5. **Review both structures.** The transform dialog identifies itself as a transformation, displays source SMILES and proposed SMILES, retains provider/model provenance and the parser-acceptance caveat, and labels the approval button “Reemplazar molécula seleccionada”. Failures and user decline never change the canvas.
6. **Use normal canvas commands in one macro.** A canvas operation wraps source deletion and proposed-graph insertion at the original component's bounding-box center inside one outer `QUndoStack` macro; it selects the new atoms after commit. Existing `DeleteSelectionCommand` removes only the selected component and its internal bonds; `_insert_molgraph_at` supplies normal add commands. Preserve supported atom/bond stereo and group metadata when the delete command restores the source on Undo. Undo restores the source molecule and removes the proposal; redo reapplies the replacement. Other components and annotations are not passed to deletion.
7. **No Clean2D integration.** The validated proposal is inserted with existing ChemIO depiction coordinates, translated to the source center. Geometry/semantic acceptance is not inferred and no M02 call, scoring gate, or chemical heuristic is added.
8. **Offline verification.** Inject fake source-SMILES export and generation functions. Test selection gating, transformed prompt content, no provider call on export failure, success review, stale-source rejection, decline/failure preservation, single-step undo/redo, and source-component replacement. Do not manually test a real provider.

## Risks / Trade-offs

- [Source graph changes during provider latency] → Compare the full source signature and selection immediately before replacement; fail closed if changed.
- [A user expects editing only a selected fragment] → Reject partial selections explicitly and state that the action requires one complete molecule; fragment-editing semantics remain out of scope.
- [Source SMILES cannot be exported] → Stop before M23/provider invocation, show a controlled failure, and keep the editor unchanged.
- [An undo macro could capture more than the molecule] → Pass only the validated component atom IDs to deletion, insert through existing canvas commands, and assert surrounding graph state through undo/redo tests.
- [A valid SMILES may be semantically wrong] → Keep the existing caveat visible and require explicit user approval; parser validation is not semantic verification.

## Migration Plan

No migration. Existing M23 generation, saved documents, provider profiles, and Clean2D evaluation remain unchanged. Removing the transformation QAction/context restores the current draw-new flow.

## Open Questions

None for the user-authorized whole-molecule scope. Structured atom/bond editing and non-destructive product insertion remain separate decisions.
