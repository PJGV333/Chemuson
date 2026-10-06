# Validation — AI/reference structure resolution

## Offline checks

- `python -m compileall -q src tests tools packaging`: passed.
- `pytest --collect-only -q`: **2019 tests collected**; no collection errors.
- `pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_recovery.py tests/test_molecular_assistant_reference_resolution.py tests/test_molecular_identity_verification.py tests/test_name2structure_service.py tests/test_name2structure_ui.py`: **138 passed**.
- `pytest -q tests/test_molecular_identity_ui_policy.py tests/test_platform_settings.py`: **18 passed**.
- `pytest -q tests/test_molecular_assistant_lifecycle.py`: **9 passed**.
- `pytest -q tests/test_molecular_assistant_transform.py`: **5 passed**.
- Focused main-window preview/mismatch/reference-fallback/undo tests: **4 passed**.
- `pytest -q tests/architecture`: **280 passed**.
- `openspec validate integrate-ai-reference-structure-resolution --strict`: valid.
- Ruff on all changed Python files with `F401,F811,F821,E722,E741`: passed.
- `git diff --check`: passed.
- The repository-wide Ruff command reports only the baseline unrelated `F401` for `math` at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- No monolithic `pytest -q` was run.

The full `tests/test_molecular_assistant_ui.py` process has a known teardown failure in the unrelated `_DescriptorWorker`: it raises `RuntimeError: wrapped C/C++ object of type _DescriptorWorker has been deleted` and Qt aborts at process exit. This failure was observed before this change as well. One assertion in that file expected the previous provenance wording; it was updated and the affected mismatch, fallback, insertion, and undo tests pass individually. The full file was not rerun after that assertion fix to avoid repeating the known Qt abort. The dedicated Molecular Assistant lifecycle suite passes.

## Optional live smoke checks

The existing local endpoint at `127.0.0.1:8080` was reachable; no server was started. With `qwen3.8-27b`, `AI + reference`, and explicit external-reference permission:

- **Ethanol:** success; AI origin verified by the existing offline `ethanol` reference (`offline-common`); ChemIO graph available.
- **Tetrandrine:** local AI request timed out; no usable reference was returned.
- **Cholesterol:** local AI request timed out; no usable reference was returned.

Bounded direct Name→Structure checks for tetrandrine and cholesterol reached the PubChem connector but returned connector errors without a graph or cache hit. No raw connector details were retained. Consequently the requested complex live fallback/cholesterol-mismatch smoke could not be demonstrated in this environment; deterministic offline mismatch and fallback tests pass.

## Scope and readiness

No Clean2D code or algorithm was changed or invoked. No model training, dataset download, new runtime dependency, or document-format change was introduced. **Merge readiness: NOT READY** until the environment can complete the complex live-reference smoke and the separate `_DescriptorWorker` teardown issue is resolved or explicitly waived. No merge was performed.
