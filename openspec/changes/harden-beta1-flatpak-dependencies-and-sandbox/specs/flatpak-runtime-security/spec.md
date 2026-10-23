## ADDED Requirements

### Requirement: Runtime image dependencies are necessary and Pillow-free

The ChemUSON runtime dependency declaration and official package manifests SHALL NOT include Pillow unless a shipped runtime path imports or calls PIL. Qt-native image annotations/rendering and RDKit's native `rdMolDraw2D` path MUST remain usable without PIL.

#### Scenario: Image annotations and native RDKit drawing do not require Pillow
- **WHEN** Pillow/PIL imports are unavailable
- **THEN** Qt `QImage` annotation clipboard/persistence behavior and the `rdMolDraw2D` native module remain available
- **AND** the runtime requirements and Flatpak manifest do not install Pillow.

### Requirement: Flatpak filesystem access is limited to ChemUSON-owned state

The Flatpak manifest MUST NOT grant the entire host home directory. Direct filesystem permissions SHALL be limited to the ChemUSON data/configuration paths currently used by the application. User-selected files SHALL continue to use the existing Qt file chooser/portal integration; this requirement does not authorize a portal or persistence refactor.

#### Scenario: Application stores settings, autosaves, and templates
- **WHEN** ChemUSON needs its existing `~/.chemuson` or `~/.config/Chemuson` state
- **THEN** the Flatpak permission covers those app-specific paths without granting arbitrary home-directory access.

#### Scenario: User selects an arbitrary document
- **WHEN** the user opens or saves a document outside ChemUSON's state directories
- **THEN** the existing file chooser/portal grants access to that selected document; broad home access is not used as a substitute.

### Requirement: Flatpak network capability is disclosed and justified

The Flatpak SHALL retain network capability only for existing configured remote Assistant and explicitly allowed external name/reference lookups. Release notes MUST disclose that using these features can send prompt, molecular, or lookup data to an external service; availability of network permission MUST NOT be described as automatic consent or provider privacy assurance.

#### Scenario: User elects to use a remote feature
- **WHEN** the user configures an Assistant endpoint or explicitly enables external reference resolution
- **THEN** Flatpak network access permits that existing feature and privacy guidance tells the user to review the destination and data sent.
