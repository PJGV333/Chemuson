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

## 5. Expanded independent reference corpus and audit baseline
- [x] 5.1 Resolve 88 representative PubChem PUG REST queries and deduplicate to 78 unique molecular structures; verify each returned connectivity SMILES against its exact input with RDKit canonicalization.
- [x] 5.2 Record each unique structure, molecular formula, PubChem CID/IUPACName/connectivity SMILES, exact pre-fix ChemName name, and baseline classification in `tests/data/chemname_iupac_reference_campaign.psv`.
- [x] 5.3 Add corpus-integrity, structure-equivalence, and stable-name tests for the 53 baseline cases assessed correct (including documented systematic/retained-name variants).
- [ ] 5.4 Reassess classifications and report reference-name exactness separately from valid systematic/retained variants after family fixes.

## 6. Stage 2 — functional groups and seniority
- [ ] 6.1 Correct branched and polycarboxylic acid parent selection/suffixes with positive and mixed-function negative cases.
- [ ] 6.2 Represent both ester components (organyl + acid-derived anion) for aliphatic and aromatic monoesters; preserve other higher-priority groups.
- [ ] 6.3 Correct branched aldehyde chain selection and multiplicative dial suffixes without degrading mixed acid/aldehyde cases.
- [ ] 6.4 Correct repeated same-class suffixes (diol, diamine, diamide, dione) and N-substituted amide locants; retain explicit mixed-priority counterexamples.
- [ ] 6.5 Add independent regression cases for chain selection, suffix/prefix priority, functional group locants, and alternate atom orders.

## 7. Stage 3 — aromatic parents and decorated substituents
- [ ] 7.1 Investigate the simple benzoic acid, benzoate ester, acetophenone, and hydroxyacetophenone `N/D` cases; fix only reference-backed general patterns.
- [ ] 7.2 Expand mono-, di-, and trisubstituted aromatic orientation/locant tests, including mixed substituents and functionalized chains.
- [ ] 7.3 Keep unsupported fused, nested, stereogenic, charged, or isotopic decorations explicitly fail-closed unless the naming path represents them fully.

## 8. Stage 4 — structural safety and determinism
- [ ] 8.1 Audit selected parent, functional groups, and rendered substituents for unaccounted heavy atoms before returning names; ensure this does not turn deliberate unsupported boundaries into blanket `N/D`.
- [ ] 8.2 Add atom-order/permutation determinism checks for symmetric and locant-tie cases.
- [ ] 8.3 Retain focused negative cases for omitted groups, duplicated groups, lost stereochemistry, charge/isotope loss, and unsupported connectivity.

## 9. Stage 5 — final coverage report
- [ ] 9.1 Classify each reference structure as correct, incorrect, unsupported, or reference pending; report overlap-aware metrics by chemical family.
- [ ] 9.2 Report discovered/corrected defects, remaining `N/D`, limited architecture impact, tests, commits, Git status, and next proposal without claiming general IUPAC conformity.
