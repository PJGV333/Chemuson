# AI Assistant Integration Specification

## ADDED Requirements

### Requirement: Runtime and development dependencies are distinct and complete
`requirements.txt` SHALL contain all runtime dependencies declared by `pyproject.toml`, including Pillow, and SHALL NOT contain test/lint/catalog tooling. `requirements-dev.txt` SHALL provide pytest, Ruff, and PyYAML. README SHALL explain a clean editable development checkout without prescribing the venv name. Open Babel SHALL remain an optional system executable.

#### Scenario: Install a new development checkout
- **WHEN** a developer follows the documented venv, editable-install, and dev-requirements steps
- **THEN** runtime and development dependencies are available without relying on the environment name `chem`.

### Requirement: Evaluator reports selected candidate provenance
For a successful validated proposal, the evaluator SHALL report selected source, selected outcome state, JSON-safe selected score when available, stable selected reason, candidate summaries, and rejected candidate summaries via the existing Clean2D summary helper. The evaluator SHALL NOT alter candidate selection or Clean2D code.

#### Scenario: Candidate is selected and others are rejected
- **WHEN** Clean2D returns selected and rejected candidates
- **THEN** JSON output identifies the selected source/state/score/reason and separates accepted summaries from rejected summaries.

### Requirement: API-key environment access is explicit
The evaluator SHALL NOT read `OPENAI_API_KEY` or any other environment variable by default. It SHALL read only the variable explicitly named by `--api-key-env`. No key may appear in report output or diagnostics.

#### Scenario: Environment contains a key without explicit flag
- **WHEN** `OPENAI_API_KEY` is set and the CLI runs without `--api-key-env`
- **THEN** provider configuration receives no API key.

### Requirement: Verification remains bounded
Automated test/diagnostic commands SHALL have a hard timeout no greater than ten minutes. A known long full suite SHALL be covered by collection and bounded shards rather than one monolithic command.

#### Scenario: Full suite exceeds the bound
- **WHEN** historical full-suite runtime exceeds ten minutes
- **THEN** verification uses collection and independent bounded test shards and never leaves a pytest command running beyond its timeout.
