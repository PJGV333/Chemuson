## ADDED Requirements

### Requirement: Molecular Assistant Has a One-Way Application Boundary

The M23 molecular-assistant package SHALL depend at runtime only on M00 Core, M01 ChemIO and the Python standard library. It SHALL NOT import M02 Clean2D, M04 ChemName, M08-M13 GUI modules, M16 name2structure, M19 bootstrap or `tools`. M00, M01 and M02 SHALL NOT import or depend on M23. AST boundary tests SHALL enforce these edges, including local and type-checking imports.

#### Scenario: M23 imports stay within Core and ChemIO
- **GIVEN** a Python source file owned by M23
- **WHEN** architecture tests analyze its imports
- **THEN** any ChemUSON dependency is limited to M00 or M01 and no forbidden module or `tools` import is present

#### Scenario: Clean2D remains usable without AI
- **GIVEN** source files owned by M02 and a catalog with M23 registered
- **WHEN** architecture tests analyze imports and dependencies
- **THEN** M02 does not import, depend on, or contact M23/provider infrastructure

#### Scenario: Importing the application contract has no graphical or parser side effects
- **GIVEN** a fresh Python process
- **WHEN** it imports `chemuson.molecular_assistant`
- **THEN** GUI, Clean2D, ChemName and RDKit modules remain unloaded and no network connection is opened

#### Scenario: The application request cannot mutate a document
- **GIVEN** the provider-neutral request contract
- **WHEN** its fields are inspected
- **THEN** it contains only the user's description and no canvas, document, selection, undo or mutation collaborator
