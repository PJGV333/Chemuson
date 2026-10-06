# Validation — GUI async worker shutdown

## Implemented scope

- Added approval-before-shutdown and deferred close handling to `ChemusonWindow`; the window tracks every running descendant `QThread` and remains alive until the last thread has stopped.
- Stopped properties/autosave timers and guarded delayed update checks and job-start entry points after shutdown approval.
- Added shutdown/result-suppression behavior to the Molecular Assistant, CompChem, template-export, and canvas-analysis owners. Name→Structure progress UI is closed before waiting; late errors/results do not update UI or documents.
- Added deterministic fake-worker tests covering Descriptor, Name→Structure, Molecular Assistant, CompChem, template export, canvas analysis, close cancellation, late-result suppression, and safe completion. No provider/model/network call is made.
- No dependencies or architecture edges changed; `architecture/modules.yml` remains unchanged.

## Passing checks

- `QT_QPA_PLATFORM=offscreen pytest -q tests/test_gui_async_worker_shutdown.py tests/test_molecular_assistant.py tests/test_molecular_assistant_ui.py tests/test_compchem3d_dock.py tests/test_geometry3d_service.py tests/test_template_controller.py tests/architecture/test_main_window_background_workers.py` — **111 passed**.
- Original crash-pair order (CompChem export then Molecular Assistant worker test) — **5/5 runs passed**; logs are `/tmp/chemuson-worker-shutdown-final-pair-{1..5}.log`.
- Reverse pair order — **5/5 runs passed**; logs are `/tmp/chemuson-worker-shutdown-final-reverse-pair-{1..5}.log`.
- `python -m compileall src tests tools packaging` — **passed**.
- `pytest --collect-only -q` — **1918 tests collected**.
- `openspec validate stabilize-gui-async-worker-shutdown --strict` — **valid**.
- `git diff --check` — **passed**.

## Baseline exceptions (not changed)

- Full `QT_QPA_PLATFORM=offscreen PYTHONFAULTHANDLER=1 pytest -q` still exits **139** at approximately 69%, while `test_controller_runs_generation_on_worker_thread_and_keeps_key_transient` processes Qt events. The C stack is the same baseline `QUndoStack` destruction / widget teardown stack recorded before implementation in `/tmp/chemuson-phase45-baseline-pytest.log`; the new minimal crash pair passes repeatedly. Full current log: `/tmp/chemuson-worker-shutdown-full-suite.log`; verbose current log: `/tmp/chemuson-worker-shutdown-full-verbose.log`.
- Before the abort, the same two baseline failures appear: `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs` and `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`. The former is independently documented in `baseline.md`; the latter passes alone and in the 111-test focused block. The prefix run confirmed both failures and otherwise passed 1055 tests (`/tmp/chemuson-worker-shutdown-prefix-to-compchem.log`). No unrelated baseline failure was altered.
- Ruff command `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` exits 1 only for the pre-existing unused `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; introduced unused imports were removed. Output: `/tmp/chemuson-worker-shutdown-final-ruff.log`.

## Architecture and external effects

- Window teardown uses the QObject descendant tree as a safety barrier. Interruption requests remain advisory; workers finish under existing bounded operations, and the GUI event loop remains responsive.
- No `QThread.terminate()`, busy wait, arbitrary sleep, external dependency, live provider, model server, or network service was used.

## Current continuation (2026-10-07)

The original CompChem→Assistant pair passed **5/5** again in both forward and reverse order. The ordered M23→transform→UI shard passed **99 tests**, and the dedicated shutdown test passed **2 tests**. A full current suite was not run because the recorded 19:26 baseline exceeds the active 10-minute cap. The complete 1946-test collection was covered by bounded shards; their precise totals, the unrelated Clean2D/stereo failures, local Qwen smoke outcomes, and the 8-minute shard split are in `../stabilize-ai-molecular-assistant-integration/validation.md`.

No second SIGSEGV trigger was reproduced. The historically recorded monolithic SIGSEGV is not claimed eliminated: no current monolithic run was made, and isolated/bounded shards cannot prove absence of an order-dependent crash outside the exercised regions.
