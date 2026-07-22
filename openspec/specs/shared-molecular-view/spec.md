# shared-molecular-view Specification

## Purpose
Define the Core-owned read-only molecular graph adapter, its supported
heterogeneous graph representations, isolation from higher-level subsystems,
and compatibility guarantees for the historical ChemName import and exception.
## Requirements
### Requirement: Core owns the molecular graph view
Core SHALL provide `MolView` and `MolecularViewNotSupported` as a read-only
adapter for the supported heterogeneous molecular graph shapes without importing
ChemName, ChemCalc or GUI modules.

#### Scenario: Core view reads a supported graph
- **WHEN** a caller constructs `MolView` with a supported MolGraph, mapping, or
  iterable graph shape
- **THEN** its atom, bond, metadata, and structural-bond queries preserve the
  established results without mutating the graph

### Requirement: Historical ChemName imports remain compatible
`chemuson.chemname.molview.MolView` SHALL be the Core `MolView` object, and
`ChemNameNotSupported` SHALL be the same exception object as
`MolecularViewNotSupported`.

#### Scenario: Historical caller catches an unsupported graph error
- **WHEN** the historical ChemName `MolView` cannot read a graph
- **THEN** the established message is raised and is catchable as
  `ChemNameNotSupported`
