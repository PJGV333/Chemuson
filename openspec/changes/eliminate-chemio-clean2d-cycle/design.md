## Context

M01 currently imports M02 for geometry, safety, coordinate application,
scaffold depiction, and block unwrap. M02 imports the M01 scoring module.
`rdkit_io.py` therefore owns both low-level import services and high-level 2D
selection. The dependency is architectural debt: five runtime exceptions and
one high-severity M01/M02 cycle.

Before:

```text
chemio.depiction_candidates -> clean2d.geometry -> clean2d.safety
clean2d.block_unwrap -> chemio.depiction_candidates
clean2d.scaffold_depiction -> chemio.depiction_candidates
chemio.rdkit_io -> clean2d.geometry/block_unwrap/scaffold_depiction
```

After:

```text
clean2d.depiction_quality -> core + clean2d internals
clean2d.block_unwrap -> clean2d.depiction_quality
clean2d.scaffold_depiction -> clean2d.depiction_quality
clean2d.imported_depiction -> depiction_quality + block_unwrap + scaffold_depiction
                               + chemio.rdkit_io/rdkit_safe
chemio -> core
```

## Goals / Non-Goals

**Goals:**

- Make M02 the single owner of imported-depiction quality, scoring, candidate
  ranking, reports, scaffold and block-unwrap orchestration.
- Preserve all observable candidate and fallback behavior.
- Make M01 depend only on M00 and retain M02 -> M01 for isolated RDKit and
  Molfile parsing services.
- Reduce active exceptions from six to the persistence TYPE_CHECKING exception
  only, with no documented cycles.

**Non-Goals:**

- Change heuristic calculations, candidate source names, metadata, scores,
  ordering, workers, timeout behavior, chemical formats, GUI behavior, or
  persistence.
- Resolve M01 -> M09 persistence typing debt.
- Move RDKit workers into Clean2D, add dependencies, or alter Core.

## Decisions

### Split imported depiction inside M02

`clean2d.depiction_quality` owns `score_imported_depiction`,
`block_donut_score`, and their private geometry helpers. It depends only on
M00, M02, and the standard library. `clean2d.imported_depiction` owns
`DepictionCandidate`, worker-result interpretation, candidate ranking, report
rows, debug output, and scaffold/block-unwrap integration. This keeps quality
pure while allowing the orchestration module to consume M01 services.

The alternative of retaining a ChemIO wrapper is rejected because every
runtime or lazy wrapper would retain M01 -> M02. Duplicating quality logic is
rejected because it risks score drift.

### Explicit ChemIO import service and scaling

`molfile_to_molgraph` will accept an optional `target_bond_length` and retain
the existing default. It performs its established `_scale_to_default` work in
ChemIO. M02 passes its requested target length and never imports the private
scaler. Extracting a Core utility is rejected because scaling is only needed
at this ChemIO format boundary in this change.

### Deliberate Clean2D surface

M02 reexports `DepictionCandidate`, `smiles_to_depiction_candidates`,
`smiles_to_molgraph_best_depiction`, and
`smiles_to_molgraph_best_depiction_with_report` from `chemuson.clean2d`.
Quality functions remain internal M02 implementation details. The old M01
paths are not preserved: M01 declares `public_api: []`, so the internal route
is an acceptable intentional break and all repository consumers migrate.

### Behavior-preserving movement

The moved function bodies, broad exception handling, sort keys, deep copies,
metadata merges, source formatting, debug output, fallback error composition,
and worker timeout arguments remain unchanged. Imports are redirected to
deliberate M01 services or local M02 modules only.

## Risks / Trade-offs

- Import-path monkeypatches can silently target an obsolete module -> migrate
  focused tests and assert public Clean2D ownership.
- M02 importing an M01 private helper would create an unstable contract -> add
  the target-length parameter to the existing public parser instead.
- Moving code can alter import initialization -> preserve lazy worker access
  where it existed and verify isolated-worker tests.
- Scoring changes would affect selected coordinates -> characterize sort order,
  scores, metadata, rejection, and non-mutation before movement.

## Migration Plan

1. Add characterization tests on the current paths.
2. Add the active OpenSpec artifacts and validate them.
3. Move quality code, then update M02 consumers.
4. Move orchestration, migrate repository imports, and remove M01 ownership.
5. Update catalog, frozen baseline, architecture tests, and documentation.
6. Run focused, structural, full-suite, Ruff, and strict OpenSpec validation.

Rollback is a revert of the ordered commits; no persisted format or external
API migration is required.

## Open Questions

None. The inspected implementation shows the parser can own optional target
scaling without moving workers or exposing a private ChemIO helper.
