## Why

Clean 2D now improves simple isolated aromatic rings, but common monosubstituted aromatics still rely on broader candidate selection. A narrow, measured extension can improve toluene/phenol/aniline-like layouts while preserving the protections around disubstituted, fused, macrocyclic, tetrandrine-like, and high-risk complex structures.

## What Changes

- Add a production candidate path for simple monosubstituted aromatic molecules in quick/publication Clean 2D only.
- Detect a single neutral 5/6-member aromatic ring with exactly one small neutral acyclic substituent and no stereo or interaction/selection-boundary risk.
- Regularize the aromatic ring and place the substituent outward from the attachment atom.
- Add focused tests for measurable improvement, chemical identity preservation, propose exclusion, and high-risk/fused/tetrandrine exclusion.

## Capabilities

### New Capabilities

### Modified Capabilities
- `clean-2d`: Clean 2D candidate selection shall support a narrow monosubstituted aromatic template candidate without weakening existing simple-aromatic, propose, complex-policy, or invariant behavior.

## Impact

- Production: likely `src/chemuson/clean2d/engine.py`.
- Tests: new monosubstituted aromatic layout coverage plus existing Clean 2D regression/baseline suites.
- Out of scope: disubstituted aromatics, fused systems, naphthalene, macrocycles, tetrandrine-like structures, high-risk complex-policy cases, propose mode, UI, public formats, RDKit/CoordGen changes, and geometry gates.
