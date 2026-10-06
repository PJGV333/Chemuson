## 1. OpenSpec and baseline
- [ ] 1.1 Verify `origin/ai/molecular-assistant-foundation` equals the expected starting HEAD and record status/baseline before implementation.
- [ ] 1.2 Define resolution modes, provenance, fallback, network policy, generation-exhaustion diagnostics, and explicit non-goals in this change.

## 2. Provider exhaustion diagnostics
- [ ] 2.1 Extend provider responses/results with allowlisted finish reason and bounded completion/reasoning token counts; never retain reasoning content.
- [ ] 2.2 Classify empty content plus `finish_reason=length` as `generation_exhausted`, preserve safe diagnostics, and prove no format-repair request is made.

## 3. Reference orchestration
- [ ] 3.1 Add typed `AI`, `AI + reference`, and `Chemical reference` modes and stable structure-origin/outcome contracts at the existing M10 controller boundary.
- [ ] 3.2 Reuse explicit name extraction and Name→Structure; revalidate reference SMILES with ChemIO and compare using isolated InChI without duplicate lookup.
- [ ] 3.3 Reconcile AI/reference matches, mismatches, AI failures, missing references, and reference-only results; preserve AI failures diagnostically.
- [ ] 3.4 Ensure open-ended requests and whole-molecule transforms never query references; no new Qt thread, browser, tool use, external dependency, or Clean2D import.

## 4. UI, settings, insertion
- [ ] 4.1 Add method selection defaulting to `AI + reference`; omit provider calls/config requirements in reference-only mode.
- [ ] 4.2 Rename/document external-reference opt-in; persist only non-secret preferences and pass only the extracted name to the resolver.
- [ ] 4.3 Show fixed provenance/source labels, both SMILES on mismatch, reference-default explicit selection, separate AI override, and controlled fallback messaging.
- [ ] 4.4 Insert either candidate through the existing undoable canvas path; prove Undo/Redo and preserve provenance through preview/insertion feedback without changing `.cmsn`.

## 5. Tests and validation
- [ ] 5.1 Add offline integration tests for the complete decision matrix, explicit-name gating, external policy, ChemIO rejection, provenance, and no model browsing.
- [ ] 5.2 Run and record Molecular Assistant, identity, Name→Structure, settings, architecture, and lifecycle tests under the 10-minute command cap.
- [ ] 5.3 Run compileall, strict OpenSpec validation, modified-file Ruff rules, and `git diff --check`; do not run monolithic pytest.
- [ ] 5.4 After offline gates only, attempt bounded live ethanol/tetrandrine/cholesterol smoke checks if the existing local Qwen endpoint and PubChem are available.
- [ ] 5.5 Commit focused changes and normally push this branch only; stop after push, with no merge/rebase/force-push.
