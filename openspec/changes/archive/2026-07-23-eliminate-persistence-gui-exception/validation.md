# Validation Report

Date: 2026-07-23

Branch: `architecture/phase5-eliminate-persistence-gui-exception`

## Baseline

Before removing the catalog row, the architecture suite reported 195 passing
tests and 3 expected failures caused by the intermediate state: production no
longer imported `ChemusonCanvas`, while `architecture/modules.yml` and two
regressions still described the historical M01-to-M09 exception. OpenSpec
validated all 18 items. Full Qt collection was already blocked by the image
dependency described below.

## Final results

| Check | Result |
| --- | --- |
| Persistence, coordination and architecture | 209 passed |
| Persistence and architecture focused set | 122 passed |
| Autosave compatibility | 21 passed |
| `python -m compileall -q src tests tools packaging` | Passed |
| `openspec validate eliminate-persistence-gui-exception --strict` | Valid |
| `openspec validate --all --strict` | 18 passed, 0 failed |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | One established F401 in `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py` |
| `git diff --check` | Passed |

The full collection discovered 1025 tests and stopped with 56 collection
errors because the Work image does not provide the system library
`libEGL.so.1`, required by PyQt6. Iteratively excluding only files that could
not import for that reason produced 1019 collected tests; execution finished
with 978 passed, 20 skipped and 21 failed. Five of those failures were obsolete
persistence fakes and were corrected to implement `rebuild_persistence_view`;
the final persistence/coordination/architecture rerun then passed 209 tests.

The remaining broad-suite failures are environmental or pre-existing and
outside this change:

- GUI/controller imports still require the missing `libEGL.so.1`.
- Isolated RDKit candidate/stereochemistry checks differ in this container
  (`rdkit_isolated` unavailable or missing expected CIP/wedge results).

No exception, skip, dependency, production fallback or baseline relaxation was
added to conceal those results.

## Architectural outcome

- M01 current and target dependencies are exactly M00.
- Every module has an empty `temporary_exceptions` list.
- Every module has an empty `circular_dependencies` list.
- `chemuson.chemio.persistence` has no M09 import in any scope.
- The former normalized M01-to-M09 exception identity remains explicitly
  rejected by the frozen-baseline regression.

## Archive gate

Task 11.1 remains intentionally open. Archive only after the current diff and
this validation report receive review.
