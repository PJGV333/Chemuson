# Baseline — 2026-10-08

- Repository: `/home/unison-pjgv/Documentos/GitHub/Chemuson`
- Branch: `ai/molecular-assistant-foundation`
- Starting local HEAD: `ab9556b5c642fc0125354e85742205cee6d23acd`
- `origin/ai/molecular-assistant-foundation`: the same SHA, verified before edits with `git ls-remote`.
- `git status --short`: no output (clean).
- `python -m compileall -q src tests tools packaging`: exit 0, no output.
- `pytest --collect-only -q`: 1968 tests collected in 0.70s, exit 0.
- Focused baseline `pytest -q tests/test_molecular_assistant.py tests/test_molecular_assistant_transform.py tests/test_molecular_identity_verification.py`: 102 passed in 13.11s.
- `pytest -q tests/test_molecular_assistant_ui.py`: 19 passed in 7.41s.
- `pytest -q tests/architecture`: 279 passed in 14.03s.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: one unrelated existing F401 at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`import math`).
- The monolithic `pytest -q` was not run: the user set a strict 10-minute per-command limit and explicitly forbids the monolithic suite. Every test command in this change remains under that limit.
