# Baseline — 2026-10-08

- Repository: `/home/unison-pjgv/Documentos/GitHub/Chemuson`
- Starting branch: `ai/molecular-assistant-foundation`
- Starting local HEAD: `bbfe663533afe4c7729eaaa13edece4e0822ab1a`
- `origin/ai/molecular-assistant-foundation`: same SHA, verified before edits using `git ls-remote origin refs/heads/ai/molecular-assistant-foundation`.
- `git status --short`: no output (clean).
- `python -m compileall -q src tests tools packaging`: exit 0, no output.
- `pytest --collect-only -q`: `1994 tests collected in 0.77s`, exit 0.
- Baseline `pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_recovery.py tests/test_molecular_identity_verification.py tests/test_name2structure_service.py tests/test_platform_settings.py`: `129 passed in 6.77s`, exit 0.
- Baseline `pytest -q tests/test_molecular_assistant_ui.py tests/test_molecular_assistant_transform.py tests/test_molecular_assistant_lifecycle.py tests/test_molecular_identity_ui_policy.py tests/test_name2structure_ui.py`: `39 passed in 21.58s`, exit 0.
- Baseline `pytest -q tests/architecture`: `279 passed in 14.49s`, exit 0.
- Baseline `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 for the known unrelated pre-existing F401 `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- No monolithic `pytest -q` was run; the user explicitly forbids it and sets a hard 10-minute command limit.
