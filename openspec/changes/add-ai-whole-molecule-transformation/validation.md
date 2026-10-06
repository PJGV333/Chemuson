# Phase 5 validation

## Changed files

- Production: `src/chemuson/gui/canvas/canvas_structure.py`, `src/chemuson/gui/command_registry.py`, `src/chemuson/gui/commands/transform_commands.py`, `src/chemuson/gui/controllers/molecular_assistant_controller.py`, `src/chemuson/gui/dialogs/molecular_assistant_dialog.py`, `src/chemuson/gui/main_window.py`, `src/chemuson/gui/main_window_ui_builder.py`, `src/chemuson/gui/shell/assembly.py`.
- Tests: `tests/test_molecular_assistant_transform.py`.
- OpenSpec: all files under `openspec/changes/add-ai-whole-molecule-transformation/`.

## Scope and architecture

- The transform uses the existing GUI/Molecular Assistant/ChemIO path (M08/M23/M01) and canvas command stack (M10). No production package or external dependency was added, and `architecture/modules.yml` was not changed.
- The ordinary draw-new flow and Phase 4.5 Clean2D evaluation remain independent; no Clean2D or chemical-acceptance behavior was changed.
- Transformation tests inject fake source-SMILES export and generator functions. No real provider, model server, or external network was used.

## Focused verification

- Focused UI, controller, canvas undo, worker lifecycle, architecture-boundary, Clean2D-geometry-contract, and Phase 4.5 tests: **138 passed**; `/tmp/chemuson-phase5-focused-tests.log`.
- M23 contract plus transformation tests, including the full-suite-order timer regression check: **72 passed**.
- `python -m compileall src tests tools packaging`: **exit 0**; `/tmp/chemuson-phase5-compileall.log`.
- `pytest --collect-only -q`: **1922 tests collected**; `/tmp/chemuson-phase5-collect.log`.
- `openspec validate add-ai-whole-molecule-transformation --strict`: **valid**.
- `git diff --check`: **passed**.
- Scoped Ruff (`ruff check src tests tools packaging --select F401,F811,F821,E722,E741`): **exit 1 only for the pre-existing F401 `math` import** at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`, identical to the recorded baseline; no new findings.

## Full-suite comparison

The final `pytest -q` run completed 1336 tests before aborting with exit 139 (SIGSEGV) near 69%. At that point 1326 tests passed, 8 were skipped, and the two failures were the already documented baseline failures:

- `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
- `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`

The abort remains in the existing Molecular Assistant/Qt teardown path: `QSignalSpy.wait()` with `QUndoStack` destruction/index-change frames. This matches the pre-Phase-5 baseline and the post-lifecycle-fix run in `baseline.md`. All four Phase 5 transformation tests ran and passed before the abort. No unrelated baseline issue was changed.

A diagnostic prefix run exposed an order-dependent test-fixture issue: the window's chemical-properties timer could recalculate deliberately decorated source metadata between snapshot and review. The test fixture now stops that timer after seeding the canvas; the final transformation tests and focused block pass with this deterministic setup.

Logs are kept outside the repository under `/tmp/chemuson-phase5-full-suite-final.log` (historical run), `/tmp/chemuson-phase5-full-prefix.log` (diagnostic identification of the fixed fixture issue), `/tmp/chemuson-phase5-focused-tests.log`, and the paths listed in `baseline.md`.

## Continuation — typed transform and Insert Variant

The request now crosses M23 as `MolecularTransformationRequest(source_smiles, instruction)` and calls `MolecularAssistant.transform()`; the generic `transform_request` hook was removed. Source export remains in the worker, and the existing strict decoder/ChemIO validator is reused. Preview offers Insert Variant, Replace Original, and Cancel. Insert Variant inserts a separate translated component with one undoable paste-style operation; its test proves Undo/Redo affects only the variant. Replace Original retains the existing atomic macro.

Current focused results: `tests/test_molecular_assistant_transform.py` **5 passed**, `tests/test_molecular_assistant_ui.py` **19 passed**, `tests/test_gui_async_worker_shutdown.py` **2 passed**, and ordered M23→transform→UI **99 passed**. Architecture suite **278 passed**. Current full collection is 1946 and was exercised as time-bounded shards; historical monolithic logs are not rerun because 19:26 exceeds the current 10-minute cap. The only Clean2D-related failure remains the independently recorded candidate test; `src/chemuson/clean2d/` remains untouched.
