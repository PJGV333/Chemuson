## Context

M23 is implemented as a provider-neutral service that returns a graph only after strict response decoding and isolated ChemIO validation. The existing GUI already has modeless/non-modal background job patterns, the command palette delegates to existing `QAction`s, and `ChemusonCanvas._insert_molgraph` inserts atoms and bonds through a `QUndoStack` macro. Phase 1 intentionally deferred observable canvas/undo/dirty-state integration to this UI phase.

The OpenAI-compatible configuration requires a base URL and model; the API key is optional. Phase 1 explicitly excluded persistent API-key settings. The UI therefore supplies explicit per-request configuration without storing credentials or request data.

## Goals / Non-Goals

**Goals:**
- Provide one small entry point in the Structure menu and existing command palette.
- Keep provider I/O and ChemIO validation off the GUI thread.
- Show a successful proposal, provenance, and a semantic caveat before asking whether to insert.
- Commit only the already validated graph using the existing single-macro undoable canvas path.
- Prove failure/decline state preservation with offline deterministic tests.

**Non-Goals:**
- Chat, history, memory, iteration, autonomous agents, tool calling, code execution, or direct AI edits to a graph.
- Persistent provider settings or API keys, provider discovery/management, more provider types, or live endpoint tests.
- Clean2D invocation, geometry optimization, edits to Clean2D, or claims that a parser-accepted graph matches the natural-language request.
- Changes to M23's request/result/provider contract or ChemIO validation route.

## Decisions

1. **Use a dedicated GUI controller with a worker thread.** The QAction calls the GUI request flow; a small M08 controller owns the per-job `QThread` and worker. The worker constructs `OpenAICompatibleProvider` and invokes `MolecularAssistant.generate`; it returns only the typed result or a controlled failure. This reuses the repository's existing Qt worker/controller pattern. A synchronous call from the window/dialog was rejected because provider I/O would stall the event loop. Network cancellation is not added to M23's public contract: closing the modeless flow marks that job as abandoned, and any late result is ignored; the provider's finite timeout still bounds the worker.

2. **Collect endpoint, model, and optional masked key for the current request only.** No credential/settings persistence is introduced. This follows Phase 1's explicit deferral of persistent keys and allows the first UI to use local or externally configured OpenAI-compatible endpoints. A saved-preferences screen was rejected as a larger settings/security scope. The optional JSON-output capability is explicit and defaults off, preserving M23's own strict decoder as authoritative.

3. **Use the same entry point and dialog flow for menu and command palette.** Add one QAction, register that QAction with the existing command registry, and add it to the Structure menu. Do not create an alternate command-palette handler or chat window.

4. **Review first, then commit using the canvas's normal insertion method.** Display the exact proposed SMILES, provider/model metadata, and the parser-acceptance disclaimer. Only an explicit Insert choice calls `_insert_molgraph(..., select_inserted=True)`. This method creates the established `Paste molecule` undo macro; no raw SMILES parser or Clean2D depiction path is added in the GUI. Automatic insertion and direct model-to-canvas mutation were rejected to preserve user control and existing undo behavior.

5. **Keep failure handling fail-closed and presentation stable.** Present only stable `status` and `reason_code`; do not show raw HTTP/exception text. On any failed result, decline, or abandoned job, do not call insertion. Tests snapshot graph chemistry, selection, undo index, and clean state around each path.

## Risks / Trade-offs

- [The user may close the dialog while blocking HTTP is still running] → The modeless request can be abandoned immediately; controller ownership keeps worker/thread references until completion, ignores late results, and relies on the finite provider deadline. Do not claim that closing aborts the HTTP request.
- [Per-request configuration is repetitive] → This is an intentional small Phase 2 boundary; a saved settings/provider-management feature needs its own later scope and security decision.
- [Parser acceptance can be mistaken for scientific correctness] → Always present the explicit semantic-correctness caveat and require approval before insertion.
- [Canvas insertion is GUI-specific and may have state regressions] → Use only the existing undo macro path and assert actual undo/redo and before/after editor snapshots in focused tests.

## Migration Plan

No data migration is required. Add the action/controller/worker and focused tests. Rollback is removal of the new action and controller integration; M23 remains usable without GUI and Clean2D remains unchanged. No persistent keys, prompts, generated payloads, or document-format changes need migration.

## Open Questions

None for this bounded phase. Endpoint selection and credentials are explicit per operation; persistent settings, additional providers, Clean2D evaluation, and structured editing remain separate phases.