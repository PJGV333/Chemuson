# Tasks — ChemName / IUPAC robustness

## 1. Scope, baseline and references
- [x] 1.1 Verify the clean campaign branch is based on `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`, without modifying the release branch or `main`.
- [x] 1.2 Capture bounded ChemName/UI tests, compileall, scoped Ruff, and exact pre-change outputs for the two reproduced cases; record the unavailable/truncated full-collection evidence accurately.
- [x] 1.3 Independently verify the systematic amide rule and the exact substituted-aryl name; record the Blue Book/PubChem distinction and retrieval dates in `corpus.md`.
- [x] 1.4 Define the smallest change boundary and fail-closed conditions in this OpenSpec before code changes.

## 2. Regression corpus and red tests
- [x] 2.1 Correct the existing `smiles_acetamide` expectation from `1-aminoethanamide` to the requested exact systematic `ethanamide`.
- [x] 2.2 Add `CC(=O)Cc1cc(N)ccc1C` with exact expected name/reference and an alternate SMILES/atom-order form to verify equivalent graph naming.
- [x] 2.3 Add negative cases proving unsupported stereochemistry/connectivity and charged/isotopic aryl decoration return `N/D` rather than a name with missing information.
- [x] 2.4 Add exact ChemName unit and packaged-smoke assertions for the corrected amide and supported substituted phenyl regression.

## 3. Minimal implementation
- [x] 3.1 Stop classifying an already-recognized amide nitrogen as an independent amine occurrence.
- [x] 3.2 Name the documented substituted phenyl group with attachment-fixed locant 1, bounded existing substituent detectors, deterministic direction selection and unambiguous parentheses.
- [x] 3.3 Add fail-closed validation for unsupported ring attachment topology and omitted stereochemical/charge/isotope/functional-group metadata.
- [x] 3.4 Keep implementation within existing `chemname` modules; do not change architecture, ChemIO/RDKit, Clean2D or unrelated naming rules.

## 4. Verification and reporting
- [x] 4.1 Run red/green focused ChemName tests and the exact acceptance subset, each with an external timeout no greater than 600 seconds.
- [x] 4.2 Run the bounded ChemName/UI shard, architecture suite, compileall, scoped Ruff, strict OpenSpec validation and `git diff --check`; record test metrics and any unrelated baseline.
- [x] 4.3 Confirm unsupported sentinels remain fail-closed and the full existing ChemName acceptance corpus has no new failures; do not run the monolithic suite or retry a timed-out command without a changed diagnostic approach.
- [x] 4.4 Update the existing roadmap campaign status with exact results and remaining boundaries.
- [x] 4.5 Commit in small, verifiable commits on `chemname/iupac-robustness` (`d7823bd` for scope/baseline and `bdfd002` for implementation/regressions); verify release branch and `main` refs are unchanged. Keep commits local; do not merge, tag, publish, or push.
