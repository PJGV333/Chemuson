# Baseline — molecular identity verification

Captured before identity-verification implementation on branch `ai/molecular-assistant-foundation` at `6d4fca961cc7695ed01d78ae8f0a1ad769388384` (fetched `origin` matched; clean tree). Shared exact command summaries and logs under `/tmp` are recorded in `../stabilize-ai-molecular-assistant-integration/baseline.md`.

- `python -m compileall src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 1922 tests collected.
- Existing M23/UI/transform and GUI-worker/CompChem focused baselines passed: respectively 86 passed/2 skipped and 6 passed.
- No identity-verification API, identity state, or semantic regression test existed at baseline.
- Full monolithic suite is intentionally not repeated: historical runtime 19:26 exceeds the operator's absolute 10-minute cap.
