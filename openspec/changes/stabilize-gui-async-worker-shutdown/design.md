## Context and baseline evidence

- Starting branch/HEAD: `ai/molecular-assistant-foundation` / `2de88ce03d786bc386bd452ea828906f5f81e325`; clean tree and identical to `origin/ai/molecular-assistant-foundation`.
- Environment: Python 3.14.7, PyQt6 6.11.0 / Qt 6.11.0, RDKit 2026.03.6, pytest 9.1.1, CachyOS x86_64 (kernel 7.2.7-1-cachyos). Full package inventory is kept only in `/tmp/chemuson-segfault-pip-freeze.txt`.
- Full `PYTHONFAULTHANDLER=1 QT_QPA_PLATFORM=offscreen pytest -vv` reproduced exit 139 near 69%. It stopped while running `tests/test_molecular_assistant_ui.py::test_controller_runs_generation_on_worker_thread_and_keeps_key_transient`; the following collected test was `test_controller_rejects_bad_endpoint_without_starting_a_job`. Detailed transcript: `/tmp/chemuson-segfault-vv.log`.
- Minimal ordered pair: `tests/test_compchem3d_dock.py::test_compchem_exports_xyz_and_inputs` followed by `tests/test_molecular_assistant_ui.py::test_controller_runs_generation_on_worker_thread_and_keeps_key_transient`. Each passes alone; this order aborts reproducibly, while reversed order passes. Pair logs: `/tmp/chemuson-segfault-minimal-pair-repeat-*.log` and `/tmp/chemuson-segfault-minimal-pair-reversed.log`.
- Faulthandler identified `_DescriptorWorker.run` blocked in `molecular_descriptors_isolated()`'s subprocess communication. A Qt event loop then destroys a window-owned QThread before completion; the equivalent main-branch probe reported `QThread: Destroyed while thread '' is still running`. A temporary probe that allowed the descriptor timer to run and explicitly waited for that QThread exited cleanly on both `origin/main` and this branch.
- This is order/lifetime dependent. The CompChem test passes in isolation. The existing Qt fixture's close/deleteLater/processEvents cycle does not coordinate pending property timers or active window-owned threads.
- `origin/main` (`4068319e9a1deee7dbf239872191a8f509c8b641`) has the same general descriptor-thread/window lifecycle pattern; the exact Molecular Assistant test is branch-only. No full suite was run on main.
- Separate baseline issue, out of scope: `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs` fails isolated on main and this branch because candidate sources omit `rdkit_isolated`; it does not cause this Qt crash. No Clean2D code will be changed.

## Goals / Non-Goals

**Goals:**
- Do not destroy any QThread while it is running.
- Complete normal document dirty/update confirmation before beginning irreversible shutdown.
- Prevent late results from touching canvases, documents, docks, dialogs, status bars, selection, or undo state.
- Include descriptor, Name→Structure, Molecular Assistant, CompChem, template SMILES export, and canvas analysis threads whose QObject parents are the window or its descendants.
- Keep the GUI event loop responsive while bounded workers finish naturally.

**Non-Goals:**
- Cancelling a blocking subprocess or HTTP operation by claiming `requestInterruption()` is sufficient.
- Changing computation, chemical semantics, worker timeouts, Clean2D, persistence, or provider contracts.
- Introducing a generic application-wide executor or refactoring unrelated workers.

## Worker ownership audit

- `ChemusonWindow`: `_DescriptorWorker`/`_descriptor_jobs` (`QThread(self)`, ChemIO isolated descriptor timeout 5s); `_NameToStructureWorker`/`_name2structure_jobs` (`QThread(self)`, connector timeout 8s).
- Window-owned controllers: `MolecularAssistantController` (`QThread(self)`, provider has its configured finite timeout) and `CompChem3DController` (`QThread(self)`, backend settings carry finite timeouts).
- Window/child-owned workers: `TemplateController._SmilesExportWorker` uses `QThread(context.parent)`; `ChemusonCanvas._CanvasAnalysisWorker` uses `QThread(self)`. Both can be descendants of the main window and need late-result suppression plus the common close barrier.
- Approved window shutdown also stops the per-document autosave timers and the pending chemical-properties timer, and guards the delayed startup update-check callback from starting during deferred close.
- `geometry3d.service._EXECUTOR` is a module-global `ThreadPoolExecutor`, not a QObject/window-owned QThread. It is not part of this teardown change; the observed `chemuson-3d_0` stack was idle and was not the fatal active descriptor QThread.
- Repository-wide `rg QThread` found no other window-owned QThread paths.

## Decisions

1. **Defer final close, rather than block the GUI.** After existing close confirmations and pending updater decisions succeed, mark shutdown, stop timers/actions that can start work, ask worker-owning controllers/canvases to suppress results and request cooperative interruption, then ignore the current close event if any descendant QThread remains active. Connect thread completion to a single recheck; when all workers have actually exited, re-enter close through a queued/safe callback and accept without repeating prompts.
2. **Use the QObject ownership tree as the final safety inventory.** Track all running `QThread` descendants of the window, including threads whose per-job dictionary was already cleared on worker completion but whose thread has not exited yet. No new worker may start once shutdown begins. Do not use a busy loop or `QThread.terminate()`.
3. **Separate result suppression from physical cancellation.** Controller/canvas callbacks discard outputs after shutdown. `requestInterruption()` is advisory; RDKit, provider, Name→Structure, and other blocking work is allowed to finish under its existing finite timeout. A result signal never authorizes a document mutation during shutdown.
4. **Preserve close cancellation semantics.** Dirty-document confirmations and pending update handling run first. If the user cancels, do not enter shutdown or abandon active work. Only an approved close enters shutdown; the final deferred close does not repeat those decisions.
5. **Do not solve this only in the pytest fixture.** Tests will still clean up widgets, but production window/timer/worker ownership must guarantee safety independently.

## Risks / Trade-offs

- [A provider or external resolver takes its full timeout] → The closing window remains alive but non-interactive until its finite operation finishes; no claim of interrupting blocking I/O.
- [A worker family is omitted from the registry] → Audit all direct and descendant QThreads and test every family; use `findChildren(QThread)` as the final barrier.
- [A worker completion callback changes UI during shutdown] → Guard result handlers and controller relays, not only dialog insertion.
- [Close is cancelled after asynchronous work exists] → Enter the shutdown barrier only after the existing confirmation/update gates pass.

## Migration Plan

No migration. If the application is closed during active jobs, it may remain pending until those bounded operations finish. Removing the coordinator restores current behavior and is the rollback path. The Atom/Bond graph, `.cmsn`, provider, and Clean2D contracts are unchanged.
