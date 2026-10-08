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

## Historial de publicación — preservar el error original

El push automático inicial de `e838766` falló con `fatal: could not read Username for 'https://github.com': No such device or address`. En ese momento, `git ls-remote` confirmaba que el remoto seguía en `915ebb27d0f4858e3fef221a8d457ba0f5f0e92f`; la rama local y la remota diferían. Este registro es histórico y no se elimina.

Posteriormente, el propietario publicó manualmente los commits `e838766` y `b4e0e66`. El fetch de esta ejecución confirmó GitHub en `b4e0e661c0062e5e7bb69d5ea0dd712ea7e8e38e`, sincronizado con la rama local. La tarea 5.5 queda marcada como completada por ese commit/push normal. Su restricción original de detenerse después del push correspondía a la fase de implementación; esta microfase posterior cuenta con autorización explícita del propietario para integrar a `main` si se satisfacen todos los gates.

## Validación actual de cierre documental

- Inicio comprobado: `origin/main`=`4068319e9a1deee7dbf239872191a8f509c8b641`; rama AI local/remota=`b4e0e661c0062e5e7bb69d5ea0dd712ea7e8e38e`.
- Comparación: `origin/main` es ancestro; 23 commits ahead, 0 behind; árbol limpio al inicio.
- `timeout 5m openspec validate integrate-ai-reference-structure-resolution --strict`: válido.
- `timeout 8m pytest -q tests/architecture`: 280 passed.
- `timeout 5m python -m compileall -q src tests tools packaging`: exit 0.
- `git diff --check`: limpio. `git diff origin/main...HEAD -- src/chemuson/clean2d/`: vacío.
- No se repitieron pruebas focalizadas Name→Structure/reference: esta microfase sólo cambia documentación. Los resultados previos (350 focalizadas, 280 arquitectura y smokes live descritos arriba) siguen siendo evidencia histórica, no resultados nuevos.
- La captura baseline registra que una única ejecución completa previa a los cambios acabó en SIGSEGV durante teardown Qt; no se repitió ni se declara resuelta. Ruff global mantiene el F401 histórico de un test Clean2D.

## Limitaciones del asistente y trabajo futuro — se mantienen

Un SMILES válido no prueba identidad; los LLM, incluido Qwen, pueden alucinar moléculas y producir JSON inválido. El backend NInfer probado no admite actualmente `response_format=json_object`; un modelo puede agotar sus tokens de razonamiento sin contenido final. Tetrandrina, vancomicina y eritromicina no están garantizadas por IA solamente. ChemIO valida la estructura, no su correspondencia semántica con el nombre; la identidad puede permanecer `UNVERIFIED`. No se afirma ningún resultado real de Qwen: el endpoint local no estaba disponible.

Siguen pendientes OAuth/cuenta GPT, un modelo químico pequeño especializado, investigación de fine-tuning, evaluación manual más amplia y futuras campañas de estrés Clean2D. La deuda Qt de `_DescriptorWorker` y el SIGSEGV histórico son independientes y no están declarados definitivamente resueltos. También persisten placeholders globales de OpenSpec ajenos y el F401 histórico del test Clean2D.

## Merge readiness actual

La revisión técnica confirmó que no hay cambios en Clean2D, M23 permanece aislado (sólo M00/M01 como dependencias propias), Clean2D no depende de M23, ChemIO mantiene importación aislada y no se alteraron contratos `.cmsn` en esta microfase. No se modifican la UI ni Command Palette en esta microfase. Pillow se conserva: `origin/main:pyproject.toml` ya lo declara en `[project].dependencies`, y el commit `6d36e3a` reconcilia `requirements.txt` con esa fuente; no es una dependencia nueva del proyecto ni se cambia aquí. M23 no introduce dependencias runtime adicionales.

El baseline monolítico de esta microfase terminó con dos marcadores `F` y aborto SIGSEGV durante teardown Qt/`QUndoStack`, exit 139 y sin resumen final. No se conocen los tests correspondientes a esos marcadores ni se atribuye causa. La deuda histórica `_DescriptorWorker`/Qt permanece abierta; este resultado no demuestra que se haya resuelto ni proporciona evidencia concreta de una regresión nueva causada por M23. Por instrucción del propietario, la deuda conocida se acepta provisionalmente para la integración; las pruebas focalizadas previas y las pruebas de arquitectura de esta ejecución pasan.

**MERGE READINESS: READY** — OpenSpec estricto, arquitectura, compileall, `git diff --check`, Clean2D intacto, rama remota sincronizada y `origin/main` ancestro de la rama AI quedan confirmados. Los límites de los modelos y la deuda Qt siguen explícitamente abiertos. La suite monolítica no se volverá a ejecutar.
