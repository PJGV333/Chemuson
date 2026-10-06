## Baseline

- [x] Confirm current branch and origin HEAD before lifecycle edits; preserve the uncommitted structured-response work.
- [x] Capture compileall, collection, focused assistant UI/transform, relevant architecture, and Ruff baseline.
- [x] Keep each command below 10 minutes; do not run monolithic pytest.

## Dialog/job lifecycle

- [x] Implement stable job-ID cleanup for all four Molecular Assistant registries and abandon pending work idempotently.
- [x] Connect job-scoped `finished` and `destroyed` callbacks without capturing/dereferencing a dead dialog.
- [x] Reorder shutdown to unregister/abandon jobs before closing live dialogs, and clear identity results too.
- [x] Prevent stale/late worker callbacks from reaching closed or destroyed dialogs; retain a successful preview only while its dialog is live.

## Regression tests

- [x] Add deterministic real-Qt tests for unstarted close, pending close, actual QObject deletion before shutdown, pending-job abandonment, completed preview close, accept/delete, five repeated open/close cycles, transform close, retry with a new ID, and window shutdown with dialog + worker alive.
- [x] Assert cleanup registries are empty and workers terminate; do not simulate deleted wrappers by raising RuntimeError.
- [x] Preserve/document the bounded claim: this resolves the concrete stale-dialog wrapper only, not all historical Qt shutdown/SIGSEGV issues.

## Validation

- [x] Run assistant UI and transform tests, relevant architecture tests, compileall, changed-file Ruff, strict OpenSpec, and diff checks within limits.
