# Validation — 2026-10-08

## Focused checks

All commands were bounded by `timeout 8m`; the monolithic suite was intentionally not run.

- `python -m compileall -q src tests tools packaging`: passed.
- Final `pytest --collect-only -q`: 1968 tests collected in 0.66s.
- `pytest -q tests/test_molecular_assistant.py`: 87 passed.
- `pytest -q tests/test_molecular_assistant_ui.py`: 19 passed.
- `pytest -q tests/test_molecular_identity_verification.py tests/test_molecular_identity_ui_policy.py`: 12 passed.
- `pytest -q tests/test_name2structure_service.py tests/test_name2structure_ui.py`: 8 passed, including offline PubChem-cache use without fetch.
- `pytest -q tests/test_platform_settings.py`: 11 passed.
- `pytest -q tests/test_molecular_assistant_transform.py`: 5 passed (controller call-path regression check).
- `pytest -q tests/architecture`: 279 passed.
- `openspec validate harden-ai-locality-and-credential-transport --strict`: valid.
- Changed-file Ruff selection (`F401,F811,F821,E722,E741`): passed.
- `git diff --check`: passed.

The strict OpenSpec validation of this change passes. `openspec validate --all --strict` reports this change valid (20 items pass) but exits nonzero for 33 unrelated existing capability specs whose Purpose sections remain placeholders; those specs are outside this change and were not edited. Existing full-tree Ruff baseline also has one unrelated F401 at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.

## UI worker shutdown observation

An early run temporarily placed two new policy-widget tests in the existing 21-test Molecular Assistant UI module. All assertions completed, then process shutdown emitted a `_DescriptorWorker` deleted-wrapper traceback and exited 134. The original 19-test module passed again at HEAD baseline and after the two new tests were isolated in `tests/test_molecular_identity_ui_policy.py`; each new test file passes independently. This task does not claim to resolve the broader SIGSEGV/worker-shutdown issue. No Clean2D code or worker-lifecycle code was changed.

## Explicit exclusions

No full monolithic pytest run, manual model evaluation, external network lookup, dataset download, model selection/training/fine-tuning, Clean2D stress campaign, or SIGSEGV campaign was performed. Historical SIGSEGV evidence remains as recorded in `docs/history/CAMPAIGNS.md` and is not a claim of definitive resolution.
