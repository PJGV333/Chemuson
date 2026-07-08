# Proposal: Stabilize Clean 2D Deterministic Baselines

## Why
Clean 2D now has a regression corpus, debug snapshots, and stable diagnostic metric semantics. The next step is to prove that the same corpus case produces the same observable record across repeated executions before any geometry improvements or stricter gates are introduced.

## What Changes
- Add a test-only deterministic baseline harness for Clean 2D corpus cases.
- Produce JSON-stable baseline records containing case metadata, result state, reason, candidate sources, metrics, and canonicalized snapshots.
- Canonicalize dictionaries, sets, lists, tuples, floats, and snapshot fields for reproducible comparison.
- Compare repeated executions using explicit metric tolerances from the existing metric contract.
- Document known observational instability rather than hiding changes in selected source, result state, or stable reason.

## Non-Goals
- No geometry improvements.
- No layout, candidate ranking, or candidate selection changes.
- No UI changes.
- No RDKit, CoordGen, or backend changes.
- No public document format changes.
- No aesthetic thresholds.
- No repairs for current known-delicate molecule families.
- No large `engine.py` refactor.

## Success Criteria
- OpenSpec validates in strict mode before implementation.
- Every corpus case can produce a JSON-serializable baseline record.
- Repeated executions of stable corpus cases compare equivalent after canonicalization.
- Float comparison uses explicit tolerances.
- Metrics remain diagnostic-only and do not recompute result state.
- Existing Clean 2D tests continue to pass.
