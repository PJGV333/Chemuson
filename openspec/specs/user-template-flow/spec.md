# user-template-flow Specification

## Purpose
TBD - created by archiving change estabilizar-flujo-plantillas-usuario. Update Purpose after archive.
## Requirements
### Requirement: Biblioteca de plantillas inicializada y persistida
Chemuson SHALL initialize a user template library when no library file exists, persist it as JSON, and keep the default user category available after loading or normalizing data.

#### Scenario: First run creates library
- **WHEN** the user template library is constructed with a path that does not exist
- **THEN** Chemuson creates a JSON library containing default categories, built-in templates, and the user category

#### Scenario: Existing partial library is normalized
- **WHEN** Chemuson loads an existing library with missing categories, duplicate categories, blank names, or repairable MOL headers
- **THEN** Chemuson normalizes the data, preserves valid templates, repairs supported MOL headers, and ensures every template category exists

#### Scenario: Platform default path is resolved
- **WHEN** Chemuson uses the default template library path
- **THEN** Chemuson stores the library under `CHEMUSON_CONFIG_HOME` when set, otherwise under the platform configuration directory

### Requirement: Plantillas se guardan desde selección o documento
Chemuson SHALL allow saving the current selection as a reusable template, or the full canvas graph when there is no selection, and SHALL reject empty content.

#### Scenario: Save selected structure
- **WHEN** the user requests saving a template while atoms are selected and provides a valid name and category
- **THEN** Chemuson stores only the selected structure as a template, persists the library, refreshes template views, and shows a success status

#### Scenario: Save full graph when nothing is selected
- **WHEN** the user requests saving a template with no selection and the canvas graph contains atoms
- **THEN** Chemuson stores the full graph as a template in the selected category

#### Scenario: Reject empty save
- **WHEN** the user requests saving a template and neither the selection nor the canvas graph contains atoms
- **THEN** Chemuson SHALL NOT create a template and SHALL notify the user that there is nothing to save

### Requirement: Categorías de plantillas gestionadas de forma consistente
Chemuson SHALL support creating, renaming, and deleting template categories while preserving valid templates and keeping menu and dock views synchronized.

#### Scenario: Create category
- **WHEN** the user creates a category with a non-blank name
- **THEN** Chemuson adds the normalized category if absent, persists the library, refreshes the template views, and shows a status message

#### Scenario: Rename category moves templates
- **WHEN** the user renames an existing category to a non-blank name
- **THEN** Chemuson updates the category list, moves templates from the old category to the new category, persists the library, and refreshes the template views

#### Scenario: Delete category preserves templates
- **WHEN** the user confirms deleting a category that contains templates
- **THEN** Chemuson moves those templates to the configured fallback user category, removes the category, persists the library, and refreshes the template views

### Requirement: Plantillas individuales gestionadas de forma consistente
Chemuson SHALL support renaming and deleting individual user templates by stable template ID.

#### Scenario: Rename template
- **WHEN** the user renames an existing template to a non-blank name
- **THEN** Chemuson updates only that template name, persists the library, refreshes the template views, and shows a status message

#### Scenario: Delete template
- **WHEN** the user confirms deleting an existing template
- **THEN** Chemuson removes only that template, persists the library, refreshes the template views, and shows a status message

#### Scenario: Unknown template ID
- **WHEN** a template operation references an unknown template ID
- **THEN** Chemuson SHALL NOT mutate any other template and SHALL surface a controlled error or no-op according to the UI action

### Requirement: Importación y exportación JSON de biblioteca
Chemuson SHALL export the current template library to JSON and import JSON libraries in merge or replace mode using normalized data.

#### Scenario: Export library
- **WHEN** the user chooses a destination file for export
- **THEN** Chemuson writes the current normalized library JSON to that file and reports success

#### Scenario: Merge imported library
- **WHEN** the user imports a JSON library in merge mode
- **THEN** Chemuson adds imported categories, skips duplicate templates by content signature, resolves template ID collisions, persists the merged library, refreshes views, and reports how many templates were added

#### Scenario: Replace library
- **WHEN** the user imports a JSON library in replace mode
- **THEN** Chemuson replaces the current library with the normalized imported library, persists it, refreshes views, and reports replacement success

#### Scenario: Invalid import file
- **WHEN** the selected import file cannot be read or parsed as a valid library source
- **THEN** Chemuson SHALL keep the current library unchanged and SHALL show an import error

### Requirement: Recuperación química de plantillas
Chemuson SHALL convert stored templates to `MolGraph` using MOL content first and SMILES as fallback when available.

#### Scenario: Convert from MOL content
- **WHEN** a stored template contains valid MOL content
- **THEN** Chemuson returns a `MolGraph` produced from that MOL content

#### Scenario: Fallback to SMILES
- **WHEN** a stored template contains MOL content that fails to convert and a non-blank SMILES fallback
- **THEN** Chemuson returns a `MolGraph` produced from the SMILES fallback

#### Scenario: No chemical content
- **WHEN** a stored template contains neither valid MOL content nor a usable SMILES fallback
- **THEN** Chemuson SHALL fail the template load with a controlled error and SHALL NOT insert an empty graph

### Requirement: Migración tolerante de plantillas MOL legadas
Chemuson SHALL migrate legacy `.mol` template files into the user library as a best-effort startup operation without blocking the application.

#### Scenario: Migrate new legacy MOL file
- **WHEN** a legacy `.mol` file exists with non-empty content and no existing template has the same name and MOL content
- **THEN** Chemuson adds it to the user category, persists the library, and reports that legacy templates were imported

#### Scenario: Skip duplicate legacy MOL file
- **WHEN** a legacy `.mol` file matches an existing template by name and MOL content
- **THEN** Chemuson SHALL NOT add a duplicate template

#### Scenario: Ignore unreadable legacy files
- **WHEN** legacy template discovery or file reading fails
- **THEN** Chemuson SHALL continue startup without raising an unhandled exception

### Requirement: Menú y dock reflejan la biblioteca
Chemuson SHALL render template categories and templates consistently in the templates menu and templates dock using stable template IDs.

#### Scenario: Refresh populated views
- **WHEN** the template views are refreshed and the library contains grouped templates
- **THEN** Chemuson displays categories and selectable template entries in both menu and dock, using the template ID as the action payload

#### Scenario: Refresh empty views
- **WHEN** the template views are refreshed and no selectable templates exist
- **THEN** Chemuson displays a disabled empty placeholder and still exposes save, category, import, and export actions in the menu where applicable

#### Scenario: Preview generation fails
- **WHEN** a template preview icon cannot be generated
- **THEN** Chemuson SHALL keep the template visible with an empty icon and SHALL NOT fail the entire refresh

### Requirement: Inserción de plantillas desde biblioteca
Chemuson SHALL insert templates selected from the menu or dock by loading the template graph from the library and delegating placement to the canvas insertion flow.

#### Scenario: Start click placement
- **WHEN** the user selects a valid template from the menu or dock
- **THEN** Chemuson enters template insertion mode with the template graph and a user-visible placement message

#### Scenario: Place immediately
- **WHEN** a valid template is requested with immediate placement
- **THEN** Chemuson inserts the template at the last scene position or visible center fallback, cancels insertion mode, and reports success

#### Scenario: Attach to atom on placement
- **WHEN** the user places a template on an existing atom and the canvas receives an attach target
- **THEN** Chemuson inserts the template graph and creates exactly one connection bond from the existing atom to the inserted template

#### Scenario: Template load fails during insertion
- **WHEN** the selected template cannot be loaded or converted to a graph
- **THEN** Chemuson SHALL NOT enter insertion mode and SHALL show a controlled load error

