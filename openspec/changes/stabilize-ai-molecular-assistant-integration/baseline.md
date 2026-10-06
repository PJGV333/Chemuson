# Baseline — AI assistant integration stabilization

Captured before production edits on `ai/molecular-assistant-foundation`, local and origin HEAD `6d4fca961cc7695ed01d78ae8f0a1ad769388384`; initial `git status --short` was empty. Fetch completed and confirmed remote HEAD exactly matched. A local recovery branch `checkpoint/ai-molecular-assistant-session-20261007` points to this SHA.

| Check | Result |
|---|---|
| `python -m compileall src tests tools packaging` | PASS, exit 0; exact output `/tmp/chemuson-session-baseline-compileall.log` |
| `pytest --collect-only -q` | PASS, 1922 collected in 0.85s; exact output `/tmp/chemuson-session-baseline-collect.log` |
| Scoped Ruff (`F401,F811,F821,E722,E741`) | Existing failure: one unused `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; exact output `/tmp/chemuson-session-baseline-ruff.log` |
| `tests/test_gui_async_worker_shutdown.py` | 2 passed in 0.93s |
| Molecular Assistant service + transform + UI modules | 86 passed, 2 skipped in 7.32s |
| CompChem dock + shutdown lifecycle | 6 passed in 1.63s |
| Prior full-suite observation | 19:26, with existing failures and a Qt SIGSEGV on an earlier baseline; this session does not run a monolithic suite under the operator's strict timeout. |

The original CompChem-export → Molecular-Assistant-worker pair passed 5/5 in forward order and 5/5 reversed (all individual pytest invocations bounded by 8 minutes). Collection around the historical ~69% crash position showed transformation/UI tests; the relevant bounded modules passed without a crash. The second trigger is not reproduced by these baseline shards and remains subject to wider shard verification.
