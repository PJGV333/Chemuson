## 1. OpenSpec and audit

- [x] 1.1 Record the baseline crash, environment, minimal order, main-branch comparison, and separate Clean2D failure.
- [x] 1.2 Validate this proposal/design/spec before production edits; strict validation also passes after implementation.
- [x] 1.3 Audit all window/descendant-owned QThreads: descriptor, Name→Structure, Molecular Assistant, CompChem, template SMILES export, and canvas analysis.

## 2. Shutdown contract

- [x] 2.1 Implement approval-before-shutdown and deferred close in `ChemusonWindow.closeEvent`; a cancelled dirty-document close must retain current jobs.
- [x] 2.2 Stop timers/new job entry points at shutdown, request interruption where appropriate, and suppress all late GUI/document results.
- [x] 2.3 Ensure every active QThread descendant completes before owner/controller/window destruction; no terminate, busy wait, or arbitrary sleeps.
- [x] 2.4 Add explicit shutdown/active-job semantics to Molecular Assistant and CompChem controllers and guard direct window/canvas/template workers.

## 3. Offline regression coverage

- [x] 3.1 Cover DescriptorWorker, Name→Structure, Molecular Assistant, CompChem, template export, and canvas analysis teardown with offline fake workers.
- [x] 3.2 Cover no late UI/document mutation, cancelled dirty close, no-worker close, worker completion, and safe QObject cleanup.
- [x] 3.3 Run the original crash pair five times in each order; run the 111-test worker/provider/geometry/UI block.

## 4. Verification and closeout

- [x] 4.1 Run focused tests, architecture tests, compileall, Ruff, strict OpenSpec, and diff check. Ruff reports only the pre-existing unused `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.
- [x] 4.2 Run the full suite once and document the unchanged baseline failure: it aborts with SIGSEGV at the same Molecular Assistant controller test and Qt `QUndoStack` destruction stack; two baseline failures are also visible before abort (Clean2D candidates and CompChem fake backend). The focused crash pair and worker-family block pass.
- [x] 4.3 Record files, gates, commit, push, final HEAD and tree state.
