# Integrate AI and reference structure resolution

## Why

A ChemIO-valid AI proposal can still be the wrong named molecule, and complex named requests can exceed a local model's output budget. ChemUSON already has a deterministic Name→Structure service, including offline names, a local PubChem cache, and an opt-in PubChem connector. The assistant should orchestrate those sources without giving the model browsing or network tools.

## Scope

- Add explicit `AI`, `AI + reference`, and `Chemical reference` resolution modes, defaulting to `AI + reference`.
- Reuse `extract_requested_molecule_name()`, `resolve_name_to_structure()`, PubChemNameConnector, its existing cache, and isolated ChemIO/InChI validation.
- Resolve references only for explicit named-molecule requests; reconcile matching and mismatching structures and allow a reference fallback after controlled AI failures.
- Carry stable provenance through preview and the normal undoable insertion route.
- Keep external name lookup opt-in and offline-by-default; send only the extracted chemical name to PubChem.
- Capture bounded finish-reason/token metadata and classify empty `finish_reason=length` responses as generation exhaustion without format repair.
- Add deterministic offline tests and an optional final live smoke check.

## Boundaries

No model browsing, tools, arbitrary URL selection, fine-tuning, dataset work, GPT OAuth, direct PubChem client, HTML scraping, Clean2D change/invocation, relaxed identity validation, or reference fallback for whole-molecule transforms. Keep the existing worker/thread lifecycle and use its worker sequentially for generation and reference resolution.

## Expected impact

The M10 Molecular Assistant controller/dialog/window, M16 identity metadata and reference reuse, M23 bounded provider/result metadata, platform preferences, OpenSpec, and focused tests. No new runtime dependency, package directory, document format, or Clean2D dependency is intended.
