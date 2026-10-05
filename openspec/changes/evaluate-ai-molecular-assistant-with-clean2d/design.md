## Context

M23 returns a graph only after strict response decoding and isolated ChemIO validation. The M02 Clean2D engine is intentionally pure: it returns candidate coordinates and diagnostics rather than mutating a graph. The repository already exposes `run_clean2d_engine` and `classify_clean2d_layout_quality`, so this phase can compose the existing contracts without changing either subsystem.

## Goals / Non-Goals

**Goals:**
- Provide an opt-in, repeatable route from a validated M23 result to M02's normal Clean2D engine.
- Report the same established visual-quality metrics on the validated graph's initial coordinates and the selected Clean2D candidate.
- Preserve source graph coordinates and chemical identity, and make the comparison JSON-safe.
- Make the tool testable offline and avoid storing prompts, credentials, or full provider responses.

**Non-Goals:**
- Adding AI imports, calls, or policy to `src/chemuson/clean2d/` or modifying its algorithms.
- Treating geometry metrics as pass/fail gates, changing Clean2D selection, or asserting that a model's chemistry is semantically correct.
- Automatic canvas insertion, UI changes, persistence, model benchmarking, provider discovery, or live endpoint tests.

## Decisions

1. **Place orchestration in `tools/`, not in either runtime module.** The evaluator is an explicit cross-domain campaign tool, not a new application capability. It may import M23 and M02; neither production package imports the tool or gains a reverse dependency. Adding a package under `src/chemuson/` is unnecessary.

2. **Fail closed before Clean2D.** Only a successful `MolecularAssistantResult` with its validated graph is eligible. A failure produces a controlled report with status/reason and no metrics; the Clean2D runner is not called. The CLI requires explicit endpoint, model, and description, so importing the tool or running the application never triggers network I/O.

3. **Use the existing engine and diagnostic surface unchanged.** Capture quality from the graph's original coordinates, call `run_clean2d_engine` with the requested bounded mode and target length, and measure the selected candidate coordinates with `classify_clean2d_layout_quality`. If the engine has no accepted candidate, `after` is `null`, while its existing state and stable reason are retained. Metrics are observational only and do not recompute or override the engine result state.

4. **Emit a small JSON report to stdout.** Include schema version, provider/model, proposed SMILES, status/reason, Clean2D state/reason, and before/after metrics. Do not include the natural-language prompt, API key, HTTP diagnostics, or raw model response. The operator may explicitly redirect stdout if they choose to retain a campaign record; the tool creates no files by default.

5. **Keep output finite and deterministic in shape.** Use the existing metric names: `quality_class`, `reason`, `crossings`, `min_nonbonded_distance`, `min_ring_degeneracy`, `length_rms_error`, `length_max_error`, `angle_rms_deviation`, `angle_max_deviation`, and `visual_score`. Non-finite/unavailable numbers serialize as `null`; no metric value is an acceptance gate. Preserve provider-returned model identity and stable Clean2D state/reason.

6. **Test through injected dependencies.** Unit tests use a fake provider and a fake Clean2D runner for call ordering, failure handling, output privacy, and serialization. One small deterministic graph exercises the real M02 engine under an explicit process timeout. No local model server, external API, or manual test is required.

## Risks / Trade-offs

- [A syntactically valid structure can still misrepresent the user's request] → The report marks ChemIO validation as syntactic/structural acceptance and makes no semantic-correctness claim.
- [A model-generated graph may have poor initial depiction coordinates] → This is expected input to the comparison; geometry findings remain diagnostic, not generation or Clean2D pass/fail gates.
- [Clean2D could mutate caller state in a future regression] → Snapshot coordinates and chemical identity around execution and test that the validated graph remains unchanged.
- [Provider failures may contain sensitive diagnostics] → Report only stable M23 status/reason codes and omit prompt, credentials, and raw response/error text.
- [The full suite exceeds the requested test-time ceiling] → Use focused tests with a hard five-minute command timeout and do not rerun the known 19-minute suite.

## Migration Plan

No migration. Add a standalone tool and focused tests. Removing the tool has no effect on the GUI, M23, M02, or saved documents.

## Open Questions

None for this diagnostic-only phase. Provider/API expansion and structured editing remain separate phases.
