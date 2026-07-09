## Why

Clean 2D now has regression reporting and review infrastructure, so it can safely take a narrow production improvement with measurable before/after evidence. Simple isolated or substituted aromatic rings are common molecules where candidate selection should prefer a safer, more regular depiction when available without touching complex/high-risk families.

## What Changes

- Improve Clean 2D production behavior for simple aromatic corpus cases by selecting a measurably better safe candidate when one exists.
- Add focused tests that prove at least one simple aromatic case improves while preserving chemical identity and protected complex behavior.
- Use baseline report before/after output to document changed corpus cases and review severity.
- Keep the change narrow, reversible, and non-format-changing; no UI, public format, macrocycle, tetrandrine, high-risk complex, or broad engine rewrite is included.

## Capabilities

### New Capabilities

### Modified Capabilities
- `clean-2d`: Clean 2D candidate selection shall allow a narrow simple-aromatic improvement while preserving identity and existing safety constraints.

## Impact

- Potential production modules: `src/chemuson/clean2d/engine.py`, `src/chemuson/clean2d/depiction_candidates.py`, `src/chemuson/clean2d/safety.py`, `src/chemuson/clean2d/quality_reporting.py`, or `src/chemuson/clean2d/complex_policy.py`, depending on the smallest safe implementation point.
- Test modules: a focused simple aromatic layout test and existing Clean 2D regression/baseline tests.
- No public document format, UI, selection metadata, chemical identity, backend integration, or complex-policy high-risk behavior is intentionally changed.
