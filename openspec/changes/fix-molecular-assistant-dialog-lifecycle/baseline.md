# Baseline — 2026-10-08

- Repository: `/home/unison-pjgv/Documentos/GitHub/Chemuson`
- Branch/HEAD/origin branch: `ai/molecular-assistant-foundation` / `ab9556b5c642fc0125354e85742205cee6d23acd`; `origin/ai/molecular-assistant-foundation` matched.
- `git status --short` (exact at lifecycle baseline):
  ```text
   M src/chemuson/gui/dialogs/molecular_assistant_dialog.py
   M src/chemuson/gui/main_window.py
   M src/chemuson/molecular_assistant/__init__.py
   M src/chemuson/molecular_assistant/limits.py
   M src/chemuson/molecular_assistant/models.py
   M src/chemuson/molecular_assistant/provider.py
   M src/chemuson/molecular_assistant/service.py
   M tests/test_molecular_assistant.py
   M tests/test_molecular_assistant_ui.py
  ?? openspec/changes/improve-ai-structured-response-recovery/
  ?? tests/test_molecular_assistant_recovery.py
  ```
  Structured-response work was present but uncommitted; no lifecycle edits existed yet.
- `python -m compileall -q src tests tools packaging`: exit 0, no output.
- `pytest --collect-only -q`: 1985 tests collected in 0.65s, exit 0.
- `pytest -q tests/test_molecular_assistant_ui.py tests/test_molecular_assistant_transform.py`: 25 passed in 15.29s.
- `pytest -q tests/architecture/test_molecular_assistant_ui_boundary.py tests/architecture/test_main_window_background_workers.py`: 5 passed in 0.13s.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: four F401 findings: existing unrelated `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`, plus unused `io`, `StructuredOutputCapability`, and `MAX_FORMAT_REPAIR_CONTENT_BYTES` imports in the in-progress `tests/test_molecular_assistant.py`.
- The monolithic `pytest -q` was not run; it is explicitly outside this task's bounded test policy. No broad Qt/SIGSEGV conclusion is drawn from the focused baseline.
