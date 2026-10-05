## Why

Phase 2 lets a user request and review a ChemIO-validated molecular proposal. There is not yet a controlled way to observe how those validated graphs behave when passed through the existing Clean2D engine. A small, explicit evaluator will connect those two capabilities for diagnostics while keeping AI out of Clean2D itself.

## What Changes

- Add an explicit command-line evaluation tool that calls M23 only when invoked with a configured endpoint, model, and description.
- Send only successful M23 `MolGraph` results to the existing, pure M02 Clean2D engine.
- Emit JSON containing provider/model provenance, proposed SMILES, the existing Clean2D result state/reason, and diagnostic quality metrics before/after.
- Add offline tests with fake providers/engine results and a bounded real-engine smoke case; no live model, network, or service is needed.

## Capabilities

### New Capabilities
- `ai-clean2d-evaluation`: Explicit, non-mutating evaluation of validated AI structure proposals with diagnostic-only Clean2D metrics.

### Modified Capabilities
- None. M23 and Clean2D public contracts and algorithms remain unchanged.

## Impact

The evaluator lives under `tools/` and can depend on M23 and M02 without creating a production module or a reverse M02 → M23 dependency. Tests exercise report contracts and non-mutation. There are no Clean2D algorithm, scoring, ranking, candidate, GUI, canvas, persistence, provider, or dependency changes.
