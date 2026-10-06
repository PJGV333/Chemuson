# Validation

- Regression reproduction uses a real Qt lifecycle: start a deliberately blocked fake job, call `dialog.deleteLater()`, process `DeferredDelete`, assert `sip.isdeleted(dialog)`, then close the window while the stale-wrapper registry would previously have called `dialog.close()`. The destroyed callback removes the stable job ID and abandons the worker first; shutdown completes without a deleted-wrapper exception.
- `pytest -q tests/test_molecular_assistant_lifecycle.py`: 9 passed in 3.10s. Cases cover unstarted close, pending abandonment/late result, direct QObject destruction before shutdown, completed preview close, accepted insertion/`WA_DeleteOnClose`, five open/close cycles, transform-context cleanup, a retry with a new job ID, and window shutdown with a live worker/dialog.
- `pytest -q tests/test_molecular_assistant_ui.py tests/test_molecular_assistant_transform.py`: 25 passed in 15.34s.
- Relevant architecture tests: 5 passed; full `tests/architecture`: 279 passed.
- Assistant/provider/recovery/identity focused tests: 113 passed.
- `python -m compileall -q src tests tools packaging`: passed; `pytest --collect-only -q`: 1994 tests collected in 0.64s.
- Strict OpenSpec validation passed for both `fix-molecular-assistant-dialog-lifecycle` and `improve-ai-structured-response-recovery`.
- Scoped Ruff on all modified Python files passed. Full repository scoped Ruff reports only the known unrelated existing unused `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- `git diff --check`: passed.
- **Bounded finding:** resolved the concrete stale `MolecularAssistantDialog` Python wrapper during `ChemusonWindow` shutdown. This does not prove that all historical Qt shutdown crashes or SIGSEGVs are resolved. No full test suite or all-worker-family shutdown campaign was run.
