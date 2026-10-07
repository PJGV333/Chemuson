# Pre-merge validation — AI/reference structure resolution

## PUG REST contract and correction

The pre-change `PubChemNameConnector` requested `IsomericSMILES,CanonicalSMILES,IUPACName` and parsed only `IsomericSMILES`/`CanonicalSMILES`. Bounded real PUG REST requests showed HTTP 200 responses for the English names with current keys `CID`, `SMILES`, `ConnectivitySMILES`, and `IUPACName`; the old connector consequently returned controlled `empty_smiles` for both tetrandrine and cholesterol.

The connector now requests `SMILES,ConnectivitySMILES,IUPACName`, chooses `SMILES` first, then accepts legacy `IsomericSMILES`, then `ConnectivitySMILES`, then legacy `CanonicalSMILES`. A missing field returns controlled `empty_smiles`. All received SMILES still pass the existing isolated ChemIO conversion before becoming a reference; ChemIO rejection leaves no graph and no cache entry. No HTTP payload was saved.

The prompt extraction is intact: `Dibuja la tetrandrina` extracts `tetrandrina`. Direct live PUG REST requests returned HTTP 404 for `tetrandrina` and Spanish `colesterol`, but HTTP 200 with usable current-property responses for `tetrandrine` and English `cholesterol`. Only these two confirmed aliases were added to the explicit `NAME_QUERY_ALIASES` table. Results retain original `query`, canonical `resolved_query`, resolved PubChem name, source, and cache status. No LLM translation or online translation was used.

## Live reference smoke

All live lookups used explicit external permission, temporary cache paths, the PUG REST API (no HTML scraping), and no retained payload files.

- `Dibuja la tetrandrina` / Chemical reference: PubChem returned a stereo-bearing SMILES; isolated ChemIO accepted it (46 atoms). Source `pubchem`, network lookup, `from_cache=false`; original query `tetrandrina`, canonical query `tetrandrine`. The Qt dialog preview appeared, insertion succeeded, Undo removed the graph, and Redo restored it.
- `Dibuja la tetrandrina` / AI + reference: an injected `generation_exhausted` AI result (no model call) followed the real external PubChem lookup; ChemIO accepted the 46-atom reference, preview showed PubChem fallback and the AI exhaustion message, and insertion/Undo/Redo passed.
- `Dibuja colesterol` / AI + reference: an injected timeout result (no model call) followed the real external PubChem lookup via the verified `colesterol` → `cholesterol` alias; ChemIO accepted the 28-atom reference, preview showed fallback, and insertion/Undo/Redo passed. No live mismatch was forced.
- The local Qwen/llama.cpp endpoint at `127.0.0.1:8080` refused connection. No server was started and no Qwen result is claimed. The fake-provider failures above exercise orchestration/UI fallback with real PubChem data; offline tests cover deterministic mismatch reconciliation.

## Offline validation

Post-change commands (each bounded below 8 minutes):

- `pytest -q tests/test_name2structure_service.py`: **13 passed**.
- `pytest -q tests/test_molecular_assistant_reference_resolution.py`: **15 passed**.
- `pytest -q tests/test_molecular_identity_verification.py`: **10 passed**.
- `pytest -q tests/test_molecular_assistant_recovery.py`: **18 passed** (including `generation_exhausted` recovery behavior).
- `pytest -q tests/test_molecular_assistant_lifecycle.py`: **9 passed**.
- `pytest -q tests/test_molecular_assistant_transform.py`: **5 passed**.
- `pytest -q tests/architecture`: **280 passed**.
- Total across these separate requested test invocations: **350 passed**. No monolithic suite was run.
- `python -m compileall -q src tests tools packaging`: passed.
- Ruff `F401,F811,F821,E722,E741` on the two modified Python files: passed.
- `timeout 5m openspec validate integrate-ai-reference-structure-resolution --strict`: valid.
- `git diff --check`: passed.
- `architecture/modules.yml` parses as YAML and the architecture test suite passes. No dependency or module-boundary change was made.
- `git diff origin/main...HEAD -- src/chemuson/clean2d`: empty.

The known `_DescriptorWorker` teardown issue is not rerun: the full assistant UI aggregation is deliberately excluded per task instructions. The dedicated assistant lifecycle suite and the bounded dialog/controller/canvas live smoke pass. This work does not claim all Qt teardown/SIGSEGV issues are resolved.

## Merge gates and delivery

Functional gates A–J pass: live PubChem structures and ChemIO validation; reference-only preview/insertion/Undo; live PubChem fallback for injected AI exhaustion/timeout; cholesterol reference fallback; deterministic offline mismatch; lifecycle; architecture; strict OpenSpec; and zero Clean2D changes. AI/model limitations and unrelated ChemUSON debt are documented separately in `docs/modules/M23-molecular-assistant.md`.

**Delivery gate K/L (clean working tree and local == remote) remains pending the required normal commit/push.** Task 5.5 is intentionally unchecked until a real commit and push have completed. No merge, rebase, or force-push has been performed.
