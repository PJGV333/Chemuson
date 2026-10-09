# BRANDING-001 — visible application name

## Purpose

Normalize user-visible application identity to the official spelling **ChemUSON** without renaming technical project or distribution identifiers.

## ADDED Requirements

### Requirement: Application presentation uses the official spelling
Application-owned visible presentation surfaces MUST use `ChemUSON`: the main-window title prefix, app bar, Qt application display name, About and Help content, user-facing dialogs/notifications, Windows installer display metadata, Linux desktop launchers, Flatpak remote/index display titles, and AppStream display metadata. The About description MUST be `ChemUSON es un editor molecular libre y de código abierto para crear, editar, visualizar y analizar estructuras químicas, diagramas y anotaciones científicas.` The main-window title MUST use `ChemUSON <version> — Editor Molecular Libre`. A descriptive title suffix and the existing version `0.3.0-beta.1` may remain.

#### Scenario: User opens the application and Help
- **WHEN** the user views the main window, app bar, About/Help dialogs, or an application notification
- **THEN** the program is presented as `ChemUSON` and the window title begins with `ChemUSON`.

#### Scenario: User installs a platform package
- **WHEN** the Windows installer or Linux desktop/AppStream metadata presents the application
- **THEN** its display name is `ChemUSON` and its existing installation/update identity remains usable.

### Requirement: Technical identities and compatibility remain unchanged
Branding changes MUST NOT rename the Python package/imports, source directories, classes/modules, executable or artifact names, repository `PJGV333/Chemuson`, Flatpak/AppStream IDs, command names, update routes/remotes, Qt `applicationName`, QSettings organization/application or preference keys, user-data paths, or `.cmsn` format/identity. Windows installer `AppId`, install directory, executable name, and artifact naming contracts MUST remain stable.

#### Scenario: Existing installation or project is reused
- **WHEN** an existing Windows installation is updated/uninstalled or an existing project/configuration is opened
- **THEN** the existing installer identity, updater path, technical identifiers, preferences, and `.cmsn` compatibility are preserved.

## Acceptance

Static and focused tests verify source strings and unchanged technical IDs. The new Windows/Linux preview packages and their installer/launcher presentation remain pending the owner's manual retest; no platform presentation is considered manually accepted yet.
