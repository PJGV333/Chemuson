# Baseline — 2026-10-08

- Repository: `/home/unison-pjgv/Documentos/GitHub/Chemuson`
- Branch: `ai/molecular-assistant-foundation`
- Starting HEAD: `0d1c97cf6849f93dc680dfd6edf0f0728764a895`
- `origin/ai/molecular-assistant-foundation`: same SHA, verified with `git ls-remote` before edits.
- `git status --short`: no output (clean).
- `python -m compileall -q src tests tools packaging`: exit 0, no output.
- `pytest --collect-only -q`: `1946 tests collected in 0.71s`, exit 0.
- Focused baseline `pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_ui.py tests/test_molecular_identity_verification.py tests/test_platform_settings.py`: `109 passed in 13.52s`.
- Architecture baseline `pytest -q tests/architecture`: `278 passed in 14.35s`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: one pre-existing unrelated F401 in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`import math`).
- The monolithic suite was not run: explicit task budget excludes it; prior campaign notes record a 19:26 runtime. No command in this task should exceed 10 minutes.
