## Context

`MolView` is a read-only adapter for heterogeneous molecular graph shapes. It
currently lives in M04 because ChemName introduced it, but M03 imports it for
formula and valence calculation. M04 also legitimately imports
`implicit_h_count` from M03, producing a catalogued cycle. Core has no upward
dependencies and already owns the structural-bond predicate required by the
adapter.

## Goals / Non-Goals

**Goals:**
- Establish one Core-owned `MolView` implementation and neutral error type.
- Preserve the exact historical ChemName import and exception-catching
  behavior.
- Remove all M03-to-M04 imports, exceptions and the resolved catalog cycle.

**Non-Goals:**
- Move `implicit_h_count` into Core or alter ChemName's M04-to-M03 dependency.
- Change graph semantics, GUI imports, chemistry calculations or error text.
- Refactor unrelated module catalog debt.

## Decisions

- Put `MolView` and `MolecularViewNotSupported` in `core/molview.py`, re-export
  them from `core.__init__`, and keep `chemname/molview.py` as a documented
  re-export. This gives Core ownership while retaining the public historical
  import without duplicate logic.
- Alias `ChemNameNotSupported` to the neutral Core error. Alias identity, not a
  subclass or wrapper, preserves existing `except ChemNameNotSupported` paths
  and the established messages.
- Change only `chemcalc.formula` and `chemcalc.valence` to import from Core.
  ChemName and GUI historical imports remain valid compatibility uses, while
  M04 retains its intentional dependency on M03.
- Update the catalog as evidence of the new graph: M00 exposes the view,
  M03 depends only on M00, M04 retains M03, and the M03/M04 cycle and two M03
  exceptions are removed. Architecture tests assert no M03 exceptions remain.

## Risks / Trade-offs

- [Exception identity changes] → Test identity and historical error catches.
- [Adapter behavior diverges during move] → Move implementation verbatim and
  run characterization tests across all supported graph shapes.
- [M03/M04 debt returns] → Assert M03 imports no M04 module and has no catalog
  exceptions.

## Migration Plan

1. Add characterization and architecture regression tests.
2. Move the implementation to Core and install the compatibility exports.
3. Migrate ChemCalc imports, update catalog and documentation.
4. Run focused, architecture, full-suite and Ruff validation; archive after
   validation. Rollback is a normal revert of the extraction commits.

## Open Questions

None.
