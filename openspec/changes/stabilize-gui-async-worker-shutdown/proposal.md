## Why

A pending chemical-properties timer can start a window-owned descriptor QThread while the window is already closing. The next Qt event loop may destroy the window and its QThread before the bounded isolated-RDKit subprocess has completed, causing `QThread: Destroyed while thread is still running` and process abort/segmentation fault. Similar ownership patterns exist for Name→Structure, Molecular Assistant, CompChem, template SMILES export, and canvas analysis jobs.

## What Changes

- Give `ChemusonWindow` an explicit shutdown state reached only after ordinary dirty-document and pending-update close decisions approve closing.
- Stop sources of new background work, suppress late results, and coordinate every active QThread in the window's QObject descendant tree.
- Defer final window close until those threads have naturally finished; interruption is advisory only and never substitutes for waiting.
- Add offline lifecycle tests for each affected worker family, cancelled close, worker-free close, and the minimal order-dependent crash reproducer.

## Capabilities

### New Capabilities
- `gui-async-worker-shutdown`: Safe, result-suppressing shutdown of all asynchronous Qt workers owned by a ChemusonWindow and its child canvases/controllers.

### Modified Capabilities
- None. Chemical algorithms, provider contracts, document serialization, and worker computation results are unchanged.

## Impact

GUI lifecycle changes in M08/M09/M10, associated tests, and `architecture/modules.yml` only if dependency edges change (not expected). No external dependency, live provider, model server, network service, `QThread.terminate()`, Clean2D algorithm, or chemical behavior change.
