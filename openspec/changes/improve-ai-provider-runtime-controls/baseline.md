# Baseline — provider runtime controls

Captured before runtime-control code changes on `ai/molecular-assistant-foundation` at `6d4fca961cc7695ed01d78ae8f0a1ad769388384` (confirmed `origin` HEAD after fetch; clean tree). The shared command record is `../stabilize-ai-molecular-assistant-integration/baseline.md`.

- `python -m compileall src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 1922 tests collected.
- Pre-change Molecular Assistant/service/transform/UI focused set: 86 passed, 2 skipped in 7.32s.
- Pre-change controller lifecycle + CompChem area: 6 passed in 1.63s.
- Repository-wide scoped Ruff has one pre-existing `F401` (`math`) in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Monolithic pytest is not run; recorded historical runtime is 19:26 and the operator's maximum is 10 minutes.
