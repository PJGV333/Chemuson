# Design

## Decisions

1. Use the integer `job_id` as the lifecycle identity. A job-scoped cleanup helper SHALL abandon work if it is still registered and remove `_molecular_assistant_dialogs`, `_molecular_assistant_results`, `_molecular_assistant_identity_results`, and `_molecular_assistant_transform_jobs` entries without dereferencing a dialog.
2. Associate `QDialog.finished` and `QObject.destroyed` callbacks with the stable ID after a controller job has started. Their closures capture the ID, not the dialog. A dialog closed before any job has no registry entry and needs no synthetic job.
3. Cleanup is idempotent so `finished`, `destroyed`, explicit abandonment, and shutdown can converge safely. A normal successful worker completion keeps its preview/result alive until the dialog finishes; a closed dialog's worker is abandoned and any later signal is ignored.
4. During shutdown, stop producers first, clean/abandon all job IDs, close only still-live dialogs, clear all transient registries, then continue existing QThread shutdown tracking. Use PyQt6's supported `sip.isdeleted` only as secondary defense if registry cleanup cannot rule out an externally destroyed wrapper; never use a blanket swallowed `RuntimeError` as the primary design.
5. Preserve existing close/dirty-document approval semantics and defer window destruction until owned workers finish. This change resolves only the concrete stale MolecularAssistantDialog wrapper and does not establish that every historical Qt/SIGSEGV shutdown issue is fixed.

## Registry semantics

A dialog may be reused for a retry after a failed attempt. Each job gets independent ID-scoped signal connections; cleanup of an older ID is harmless and cannot remove the newer ID's state. Successful preview state remains registered until acceptance/decline/close. A failure can release that job's registry after presenting the failure, so a later retry starts without stale per-job entries.

## Verification

Exercise real `close()`/`deleteLater()`/Qt event processing, not synthetic RuntimeErrors. Test all seven requested cases with bounded fake workers, assert every registry is empty at the appropriate lifecycle point, and assert the controller has no active jobs after pending work drains. Re-run assistant UI/transform and relevant architecture tests under the 10-minute command cap.
