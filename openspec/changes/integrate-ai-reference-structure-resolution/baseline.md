# Pre-merge baseline — AI/reference structure resolution

- Repository: `/home/unison-pjgv/Documentos/GitHub/Chemuson`
- Branch: `ai/molecular-assistant-foundation`
- Expected starting HEAD: `915ebb27d0f4858e3fef221a8d457ba0f5f0e92f`
- Observed local HEAD: `915ebb27d0f4858e3fef221a8d457ba0f5f0e92f`
- Observed `origin/ai/molecular-assistant-foundation`: same SHA after `git fetch origin --prune`.
- `git status --short`: no output (clean).
- `git rev-list --left-right --count main...HEAD`: `0 21` (21 ahead of main, 0 behind).
- A local recovery ref `backup/ai-molecular-assistant-premerge-915ebb2` points to the starting HEAD.

## Baseline commands

- `python -m compileall -q src tests tools packaging`: exit 0, no output.
- `pytest --collect-only -q`: `2019 tests collected in 1.12s`, exit 0.
- `pytest -q tests/test_name2structure_service.py`: `5 passed in 0.12s`.
- `pytest -q tests/test_molecular_assistant_reference_resolution.py`: `15 passed in 7.84s`.
- `pytest -q tests/test_molecular_identity_verification.py`: `10 passed in 6.07s`.
- `pytest -q tests/test_molecular_assistant_recovery.py`: `18 passed in 0.26s`.
- `pytest -q tests/test_molecular_assistant_lifecycle.py`: `9 passed in 4.45s`.
- `pytest -q tests/test_molecular_assistant_transform.py`: `5 passed in 8.20s`.
- `pytest -q tests/architecture`: `280 passed in 15.99s`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 for the known unrelated `F401` `math` in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` only.
- No monolithic `pytest -q` was run (explicitly prohibited for this task).

## Initial live PUG REST observation (no payload retained)

- PUG REST requests for `tetrandrine` and `cholesterol` using `SMILES,ConnectivitySMILES,IUPACName` returned HTTP 200 with those current property names; both `SMILES` converted through isolated ChemIO (46 and 28 atoms respectively).
- Requests for the Spanish prompt names `tetrandrina` and `colesterol` returned HTTP 404; `extract_requested_molecule_name("Dibuja la tetrandrina")` correctly extracts `tetrandrina`.
- The existing connector's legacy request/parser returned controlled `empty_smiles` for both `tetrandrine` and `cholesterol`, despite the current PUG REST response containing `SMILES` and `ConnectivitySMILES`.
- No local Qwen/llama.cpp endpoint was listening at `127.0.0.1:8080` (connection refused). No model was started.
- Only concise status/property-name/conversion observations are recorded here; HTTP response bodies were not saved.
