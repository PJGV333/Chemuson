## Why

ChemIO currently owns imported-depiction quality, candidate ranking, and layout
orchestration even though those are Clean2D responsibilities. That placement
creates five M01-to-M02 runtime exceptions and a high-severity M01/M02 cycle,
while `rdkit_io.py` mixes low-level format I/O with 2D depiction selection.

## What Changes

- Move imported-depiction scoring and donut-quality evaluation to Clean2D.
- Move SMILES depiction candidate orchestration, ranking, reports, scaffold,
  and block-unwrap integration to Clean2D.
- Keep ChemIO responsible for format parsing, isolated RDKit workers, and
  low-level SMILES/Molfile services.
- Remove `chemio.depiction_candidates` without a compatibility shim and
  migrate all repository consumers to the deliberate Clean2D API.
- Preserve candidate sources, score values, metadata, ordering, fallbacks,
  errors, coordinate scaling, worker timeouts, and stereo behavior.
- Remove all five runtime M01-to-M02 exceptions and the M01/M02 cycle. Retain
  only the existing M01-to-M09 TYPE_CHECKING persistence exception.

## Capabilities

### New Capabilities

- None.

### Modified Capabilities

- `architecture-boundaries`: Prohibit every M01-to-M02 import and remove the
  documented M01/M02 cycle while allowing M02 to consume ChemIO import services.
- `module-catalog`: Record M01 as Core-only, M02 as owning imported depiction,
  and reduce the temporary-exception baseline to the persistence typing debt.
- `clean-2d`: Establish Clean2D ownership of imported-depiction quality,
  candidate selection, ranking, reports, and its public import path.

## Impact

Affected modules are `chemio` (M01), `clean2d` (M02), focused depiction tests,
the module catalog and architecture tests, and their module documentation.
This intentionally breaks the undeclared internal import paths
`chemuson.chemio.depiction_candidates` and depiction-selection functions in
`chemio.rdkit_io`; M01 has no declared public API. No GUI behavior, persistence
format, RDKit worker ownership, chemical heuristic, or external dependency
changes are in scope.
