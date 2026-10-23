## MODIFIED Requirements

### Requirement: Runtime and development dependencies are distinct and complete

`requirements.txt` SHALL contain all runtime dependencies declared by `pyproject.toml` and SHALL NOT contain test/lint/catalog tooling. Pillow SHALL NOT be declared or packaged unless a shipped runtime path imports or calls PIL; Qt-native image annotations and RDKit's native `rdMolDraw2D` path SHALL work without PIL. `requirements-dev.txt` SHALL provide pytest, Ruff, and PyYAML. README SHALL explain a clean editable development checkout without prescribing the venv name. Open Babel SHALL remain an optional system executable.

#### Scenario: Install a new development checkout
- **WHEN** a developer follows the documented venv, editable-install, and dev-requirements steps
- **THEN** runtime and development dependencies are available without relying on the environment name `chem`
- **AND** installing the documented runtime dependencies does not install unused Pillow.
