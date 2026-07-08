## Why

Clean 2D quality diagnostics are currently scattered across multiple modules (engine.py, safety.py, depiction_candidates.py, Clean2DController). Each module uses its own diagnostic format and vocabulary, causing duplicated logic, fragile tests, and difficulty in understanding overall quality. A narrow reporting contract with lightweight helpers will make tests more stable without changing layout behavior or forcing a cross-cutting migration.

## What Changes
- Introduce a single `clean-2d-quality-reporting` capability that defines the reporting contract.
- Map all current visual quality metrics (crossings, collisions, ring distortion, etc.) to a stable vocabulary of reasons already defined in document-clean2d-decision-contract.
- Unify how candidates are reported: source module, score/quality value, final state, and rejection reason.
- Preserve internal metadata as diagnostic context while exposing only the stable contract for tests.
- Add minimal adapters/wrappers only where a clear diagnostic surface already exists.
- Record conceptual duplication between safety checks, ranking logic, depiction candidate scoring, and controller diagnostics as follow-up work when it requires refactor.

## Capabilities
### New Capabilities
- `clean-2d-quality-reporting`: Provides a stable reporting contract for Clean 2D quality diagnostics, including state, reason, reporting score, source, candidate information, and internal metadata. Tests will validate against this contract rather than implementation details.

### Modified Capabilities
- None (no existing spec changes required).

## Impact
- Code touching Clean 2D diagnostics: contract helpers, tests, and minimal adapters near existing diagnostic surfaces in engine.py, safety.py, depiction_candidates.py, or Clean2DController as needed.
- No change to layout algorithms, geometry, RDKit/CoordGen, MolGraph, canvas, or ranking/selection behavior.
- No global refactor of engine.py.
- Minor UI updates only if needed to pass through stable diagnostics.
