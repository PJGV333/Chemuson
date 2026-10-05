# Phase 4 validation record

## Baseline

- Branch `ai/molecular-assistant-foundation`; tree clean at `5f8f5386ab6dea4b084c7e00597ee34d533728c6`, matching `origin`.
- Baseline compileall passed; 1,910 tests collected; focused provider/UI/boundary tests `85 passed, 2 skipped` in 4.85s.
- Repository Ruff retains the known unrelated baseline F401 `math` import in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Full pytest was not run. The captured full suite takes 19:26 and contains one unrelated pre-existing CompChem failure, exceeding the requested 5–10 minute maximum.

## Focused verification

- `timeout 300s pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_ui.py tests/architecture/test_molecular_assistant_boundary.py tests/architecture/test_molecular_assistant_ui_boundary.py`: **91 passed, 2 skipped** in 4.80s.
- `timeout 300s pytest -q tests/architecture/test_molecular_assistant_boundary.py tests/architecture/test_molecular_assistant_ui_boundary.py tests/architecture/test_import_boundaries.py tests/architecture/test_module_catalog.py tests/architecture/test_public_api_exists.py`: **125 passed** in 7.25s.
- `python -m compileall src tests tools packaging`: **passed**.
- `pytest --collect-only -q`: **1,916 tests collected** in 0.79s.
- Scoped Ruff on changed Python files: **passed**.
- `openspec validate add-ai-molecular-assistant-provider-profiles --strict`: **valid**.
- `openspec validate --all --strict`: **47 passed, 0 failed**.
- `git diff --check`: **passed**.

## Architectural outcome and limits

The profile catalog is part of M23 and adds no dependency. M10 remains the only GUI layer importing M23; profile metadata is injected into the dialog via the controller to preserve that boundary. Profiles for OpenAI, LM Studio, and llama.cpp use the existing strict Chat Completions adapter; model IDs are user-supplied. Fake-transport tests establish request-contract compatibility only. No live service, external API, manual test, server process, new external dependency, native vendor protocol, or persisted credential was used.

An initial architecture test caught an M23 import in the dialog; this was corrected by injecting profile metadata through the M10 controller, and the final boundary suite passes.
