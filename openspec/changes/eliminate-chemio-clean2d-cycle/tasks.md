## 0. Baseline

- [x] 0.1 Capture architecture, focused depiction, collection, full-suite, and Ruff baselines. Files: `/tmp/chemuson-phase4-baseline-ruff.json`. Acceptance: baseline is green with only the historical F401.

## 1. Reconnaissance

- [x] 1.1 Inventory M01/M02 symbols, imports, consumers, monkeypatch paths, scaling, workers, and catalog state. Files: `/tmp/chemuson-phase4-cycle-recon.md`. Acceptance: report records owners, risks, movement plan, and confirms M01 public API is empty.

## 2. Characterization

- [x] 2.1 Characterize imported depiction candidate fields, ranking, scores, metadata, errors, fallbacks, scaling, non-mutation, and determinism. Files: focused depiction tests. Acceptance: tests cover behavior without requiring a live RDKit worker when fakes are available.

## 3. Depiction Quality Ownership

- [x] 3.1 Move the single implementation of imported depiction scoring and donut quality to `clean2d.depiction_quality`. Files: `depiction_quality.py`, `block_unwrap.py`, `scaffold_depiction.py`. Acceptance: M02 consumers use the new owner and focused tests preserve scores and metadata.

## 4. Imported Depiction Orchestration

- [x] 4.1 Move candidate construction, worker-result handling, ranking, reports, and derived candidate orchestration to `clean2d.imported_depiction`. Files: `imported_depiction.py`, `rdkit_io.py`. Acceptance: M02 calls only explicit M01 format/worker services and preserves all fallbacks.

## 5. Call-Site Migration

- [x] 5.1 Migrate all repository consumers and monkeypatches to deliberate Clean2D import paths. Files: focused tests and actual consumers. Acceptance: no source, test, or tool imports `chemio.depiction_candidates`.

## 6. ChemIO Cleanup

- [x] 6.1 Remove depiction orchestration and all Clean2D imports from ChemIO, then delete `chemio/depiction_candidates.py`. Files: `rdkit_io.py`, deleted module. Acceptance: AST inventory reports zero M01-to-M02 imports and no compatibility shim remains.

## 7. Catalog and Exception Baseline

- [x] 7.1 Update the module catalog and architecture tests for cycle-free M01/M02 ownership and the one-row exception baseline. Files: `architecture/modules.yml`, architecture tests. Acceptance: M01 is M00-only, M02 is M00/M01, no cycles exist, and reintroduced historical identities fail.

## 8. Documentation

- [x] 8.1 Document ChemIO low-level ownership and Clean2D imported-depiction ownership. Files: architecture and M01/M02 module docs. Acceptance: documentation records 6-to-1 exceptions, zero runtime debt, zero cycles, and remaining persistence typing debt.

## 9. Validation

- [x] 9.1 Run compile, focused tests, architecture tests, structural AST inventory, full suite, Ruff, type-checker discovery, and strict OpenSpec validation. Files: all changed files. Acceptance: all required commands complete successfully except the one historical global Ruff F401.

## 10. Archive

- [ ] 10.1 Archive `eliminate-chemio-clean2d-cycle` only after review approval. Files: OpenSpec archive. Acceptance: intentionally pending for this phase.
