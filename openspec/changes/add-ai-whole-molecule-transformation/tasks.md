## 1. OpenSpec and baseline

- [x] 1.1 Capture clean local baseline commit, tree, compileall, test collection, full-suite baseline failure, and Ruff finding before implementation.
- [x] 1.2 Strictly validate this proposal, design, and capability specs before production edits.
- [x] 1.3 Confirm architecture boundary: reuse existing M01/M08/M10/M23 edges and introduce no production package/dependency.

## 2. Whole-molecule request flow

- [x] 2.1 Add one Structure-menu/command-palette action for transforming the complete selected molecule.
- [x] 2.2 Reject empty, partial, multi-component, or otherwise ambiguous selections without starting a worker.
- [x] 2.3 Add an optional per-request transformation hook to the existing Molecular Assistant controller; export source SMILES and compose the request inside its worker.
- [x] 2.4 Snapshot source graph/selection/document context and revalidate it before replacement; suppress or abandon safely when the dialog closes.
- [x] 2.5 Extend the review dialog to show source and proposed SMILES and require explicit “Replace” approval.
- [x] 2.6 Implement selected-component deletion plus proposal insertion as one undoable canvas macro, centered on the source component; preserve other graph components.
- [x] 2.7 Preserve all supported source atom/bond stereo and group metadata when Undo restores the replaced component.

## 3. Deterministic offline regression tests

- [x] 3.1 Cover discoverability, complete-component selection gate, and transformed prompt/source handling with fake SMILES export/generator.
- [x] 3.2 Cover provider/export failures, user decline, stale/missing source document, and unchanged editor state on each no-op path.
- [x] 3.3 Cover approved replacement and exact Undo/Redo as one stack step, preserving unrelated molecules and canvas state.
- [x] 3.4 Verify standard generate-new flow, M23 contracts, architecture boundaries, and Clean2D independence remain unchanged.

## 4. Verification and closeout

- [x] 4.1 Run focused UI/controller/canvas/architecture tests, compileall, collection, scoped Ruff, strict OpenSpec, and diff checks.
- [x] 4.2 Run a local offline worker/provider test block; never invoke a real model/provider or external network.
- [x] 4.3 Run the full suite; compare any abort/failures to the captured baseline, investigate new failures, and do not alter independent baseline issues.
- [x] 4.4 Record files, architecture decision, validation outcomes, commit, push status, final HEAD, and tree state.
