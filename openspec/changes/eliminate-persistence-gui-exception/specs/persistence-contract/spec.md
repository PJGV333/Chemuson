## ADDED Requirements

### Requirement: GUI-Neutral Persistence Document

The persistence document SHALL provide `model`, `get_persistence_data`,
`load_persistence_data`, and `rebuild_persistence_view` through a structural
contract owned by ChemIO.

#### Scenario: Structural fake is accepted
- **WHEN** a pure fake provides those members without inheriting the Protocol
- **THEN** `PersistenceManager.save_to_dict` and `load_from_dict` operate on it

### Requirement: CMSN Save Compatibility

`save_to_dict` and `save_to_file` SHALL preserve the CMSN payload, application,
version, model and canvas structures, and atomic temporary-file replacement.

#### Scenario: File is saved
- **WHEN** a document is saved to a CMSN file
- **THEN** JSON equals its dictionary payload and `os.replace` finalizes it

### Requirement: CMSN Load Compatibility

`load_from_dict` and `load_from_file` SHALL preserve legacy defaults, errors,
model restoration, next IDs and canvas payload delivery.

#### Scenario: Valid payload is restored
- **WHEN** a CMSN payload is loaded
- **THEN** the model and canvas payload are restored before one visual rebuild

### Requirement: Public Rebuild Hook

The GUI canvas SHALL expose `rebuild_persistence_view` and SHALL delegate it to
the existing visual model reconstruction implementation.

#### Scenario: Rebuild is requested
- **WHEN** persistence calls the public rebuild hook
- **THEN** the existing private reconstruction executes exactly once

### Requirement: Import Isolation

Importing `chemuson.chemio.persistence` SHALL not load `chemuson.gui` or PyQt6.

#### Scenario: Fresh process imports persistence
- **WHEN** a subprocess imports the persistence module
- **THEN** no GUI or PyQt6 module is loaded
