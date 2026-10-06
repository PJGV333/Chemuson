# Fix Molecular Assistant dialog lifecycle

## Why

`MolecularAssistantDialog` uses `WA_DeleteOnClose`, while `ChemusonWindow` keeps Python wrappers in per-job registries. If Qt destroys the C++ dialog first, a later `dialog.close()` or `dialog.job_id` access can raise `RuntimeError: wrapped C/C++ object ... has been deleted` during shutdown. A concrete user report reproduced this twice.

## Scope

- Tie transient dialog/result/identity/transform registry cleanup to stable job IDs, not a dialog lookup during teardown.
- Connect a job's `finished` and `destroyed` signals to idempotent cleanup that abandons late work and removes all associated references.
- Close only live dialogs after unregistering/abandoning jobs; suppress results during shutdown and include the identity-results registry in cleanup.
- Add deterministic Qt lifecycle regressions for unstarted, pending, previewed, accepted, repeated, transform, and shutdown cases.

## Out of scope

No broad Qt worker refactor, claim that historical SIGSEGVs are solved, HTTP cancellation, Clean2D changes, or changes to Molecular Assistant chemistry behavior. This is an incidental, narrow correction discovered during the structured-response task.

## Likely impact

`src/chemuson/gui/main_window.py`, `src/chemuson/gui/shell/assembly.py` only if registry typing/ownership requires it, the existing M08 Qt UI tests, and a small OpenSpec delta. No new dependency or cross-module architecture edge is intended.
