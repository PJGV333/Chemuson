# Design

## Architecture and trust boundaries

Keep orchestration in the existing M10 `MolecularAssistantController` worker. It already runs generation, isolated ChemIO work, and identity verification off the GUI thread, and M10 already has a catalogued dependency on M16 Name→Structure. Introduce typed resolution-mode/origin/outcome contracts in that controller boundary; do not make M23 import Name→Structure and do not add a parallel PubChem implementation.

For `AI + reference`, run the existing AI generation first, then—only when `extract_requested_molecule_name()` accepts the entire prompt—call the injected/default Name→Structure resolver with that exact name and explicit `allow_network`. Sequential execution reuses the existing QThread, gives the provider and resolver their own finite deadlines, avoids Qt lifecycle risk, and leaves the GUI responsive. `AI` never invokes the resolver. `Chemical reference` never constructs/calls the provider and rejects requests without an explicit extracted name.

Use `NameToStructureResult` as the reference contract. Revalidate its returned SMILES through isolated ChemIO before it can become a candidate; compare identity via the existing `verify_molecular_identity()` and isolated InChI path, injecting the already-resolved result to avoid a second lookup. The model receives no resolver result, PubChem credentials, browsing/tool contract, or document contents. PubChem receives only the extracted name through the existing connector. Offline mode still allows StaticNameConnector and the existing local PubChem cache; only a cache miss can perform a request when the persisted opt-in is true.

## Reconciliation and provenance

Use a closed origin enum: `ai`, `reference`, `ai_verified_by_reference`, and `ai_mismatch_reference`. The typed orchestration outcome retains the AI result (including stable failure diagnostics), validated reference result, identity result, and candidates. It does not persist prompts, model payloads, reasoning, or credentials.

- AI success + same isolated InChI: offer the unchanged AI graph as the normal candidate, mark verified, and display the PubChem/offline reference provenance.
- AI success + different InChI: do not choose or mutate either graph automatically. Show both SMILES, recommend/reference-default the reference action, and expose a separate explicit AI-override action retaining the existing confirmation guard.
- AI failure + valid reference: preserve the AI status/reason diagnostically, present a reference-origin preview, and offer reference insertion or cancellation.
- AI success + no usable reference: keep the AI graph available with identity explicitly unverified.
- AI failure + no usable reference: return the ordinary controlled failure.
- Reference-only: skip provider setup and model I/O; require a named request and show reference provenance.

Both AI and reference insertion select the intended graph and call the same `_insert_molgraph(..., select_inserted=True)` path. Provenance remains explicit in the outcome, preview, and insertion confirmation/status; no new `.cmsn` serialization fields are introduced. Transform requests force/require AI-only mode and never enter reference resolution.

## Generation exhaustion

Extend `ProviderResponse` with an allowlisted finish reason and bounded non-negative integer completion/reasoning token counts. Extract only `finish_reason`, `usage.completion_tokens`, and supported numeric reasoning-token usage fields. Never read, retain, log, display, or persist `reasoning_content`. If content is empty and finish reason is `length`, return stable `generation_exhausted` diagnostics and skip format repair; reference fallback may then run for a named request. All other response parsing, schema checks, one-shot `invalid_json` repair, ChemIO validation, and transform behavior remain strict.

## UI and preferences

Expose a method selector defaulting to `AI + reference` with `AI` and `Chemical reference` alternatives. Hide/disable provider configuration for reference-only mode. Rename the existing external permission to clearly mean chemical reference lookup via PubChem and explain that only the extracted name is sent. Persist only this non-secret option and the selected non-secret method; retain offline/cache resolution when external access is off. Use fixed UI messages and allowlisted source labels; do not display raw connector exceptions or arbitrary model strings as control values.

## Testing and risks

Use fake providers and fake reference resolvers; never access PubChem from automated tests. Cover same identity with distinct SMILES, mismatch choices, fallback for malformed/timeout/exhausted AI, absent references, network-policy propagation, name-intent gating, transform isolation, reference ChemIO rejection, provenance, undo/redo, and lifecycle. Keep Clean2D imports absent and ensure no browser/tool API exists. Run the lifecycle suite and bounded focused/architecture/OpenSpec/Ruff/compile checks. Only after offline verification, optionally perform local-provider ethanol, live PubChem tetrandrine, and cholesterol mismatch smoke checks.

Risks: the sequential reference lookup starts after AI completion, so fallback latency includes the configured model timeout; this is accepted for phase one to avoid a second Qt worker. Do not claim network cancellation aborts in-flight HTTP. A live result is not required for offline tests, but merge readiness requires the requested manual tetrandrine and cholesterol smokes when connectivity/server access permit.
