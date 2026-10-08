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
- For the original pre-implementation baseline at SHA `915ebb27d0f4858e3fef221a8d457ba0f5f0e92f`, no monolithic `pytest -q` was run, as specified at that time. This later documentation microphase ran its required pre-change baseline once; see below.

## Initial live PUG REST observation (no payload retained)

- PUG REST requests for `tetrandrine` and `cholesterol` using `SMILES,ConnectivitySMILES,IUPACName` returned HTTP 200 with those current property names; both `SMILES` converted through isolated ChemIO (46 and 28 atoms respectively).
- Requests for the Spanish prompt names `tetrandrina` and `colesterol` returned HTTP 404; `extract_requested_molecule_name("Dibuja la tetrandrina")` correctly extracts `tetrandrina`.
- The existing connector's legacy request/parser returned controlled `empty_smiles` for both `tetrandrine` and `cholesterol`, despite the current PUG REST response containing `SMILES` and `ConnectivitySMILES`.
- No local Qwen/llama.cpp endpoint was listening at `127.0.0.1:8080` (connection refused). No model was started.
- Only concise status/property-name/conversion observations are recorded here; HTTP response bodies were not saved.

## Microfase documental — baseline al 2026-10-07

- `git fetch origin --prune`: completado.
- `git status --short --branch`: `## ai/molecular-assistant-foundation...origin/ai/molecular-assistant-foundation` (limpio).
- `HEAD` y `origin/ai/molecular-assistant-foundation`: `b4e0e661c0062e5e7bb69d5ea0dd712ea7e8e38e`.
- `origin/main`: `4068319e9a1deee7dbf239872191a8f509c8b641`.
- `git rev-list --left-right --count origin/main...origin/ai/molecular-assistant-foundation`: `0 23`; `origin/main` es ancestro de la rama AI.
- `python -m compileall src tests tools packaging`: exit 0.
- `pytest --collect-only -q`: `2027 tests collected in 0.88s`, exit 0.
- Baseline completo ejecutado una sola vez antes de cambios: `timeout 10m pytest -q`. Mostró dos marcadores `F` antes de abortar alrededor del 63% por SIGSEGV durante teardown Qt (`QUndoStack`/`QWidget`, con worker `chemuson-3d_0`); terminó con exit 139 y sin resumen final. No se repitió. La deuda histórica de teardown/SIGSEGV permanece sin declarar resuelta.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó.
- Validación focal posterior del árbol de campaña: OpenSpec estricto válido; `pytest -q tests/architecture` 280 passed; compileall quiet exit 0; `git diff --check` limpio.
- `git diff origin/main...HEAD -- src/chemuson/clean2d/`: vacío.
- Decisión del propietario: conservar `Pillow`. `origin/main:pyproject.toml` ya declara el paquete en `[project].dependencies`; `6d36e3a` reconcilió `requirements.txt` con esa declaración. La ausencia de imports PIL directos no invalida la dependencia; no se modifica ninguna dependencia en esta microfase.
