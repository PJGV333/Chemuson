## Why

ChemCalc imports the molecular graph adapter from ChemName, creating an M03 to
M04 runtime dependency that closes a documented M03/M04 cycle. The adapter is
generic Core infrastructure and can be owned by Core without changing the
historical ChemName import contract.

## What Changes

- Move the single `MolView` implementation and its neutral unsupported-graph
  error into Core.
- Preserve `chemuson.chemname.molview.MolView` and
  `ChemNameNotSupported` as compatibility exports of the Core implementation.
- Change ChemCalc to depend only on Core, removing its two M03-to-M04 runtime
  exceptions and the M03/M04 catalog cycle.
- Document the shared view public API, Core ownership, and completed debt
  reduction in the architecture catalog and module documentation.

## Capabilities

### New Capabilities
- `shared-molecular-view`: Core-owned read-only molecular graph adapter with a
  stable historical ChemName compatibility import.

### Modified Capabilities
- `module-catalog`: Record the completed M03/M04 decoupling and the Core public
  API without changing the six M01 exceptions.
- `architecture-boundaries`: Prevent reintroduction of M03 temporary exceptions.

## Impact

Affected code is limited to `core`, `chemcalc`, and the ChemName compatibility
module, plus architecture tests, catalog and documentation. GUI and ChemName
callers retain their historical imports; ChemName continues to depend on
ChemCalc for `implicit_h_count`.
