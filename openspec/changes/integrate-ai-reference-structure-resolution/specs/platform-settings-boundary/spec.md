## ADDED Requirements

### Requirement: Reference network permission is a non-secret preference
The platform settings boundary SHALL persist the user's external chemical-reference permission and selected resolution method as non-secret preferences, default to external access disabled and AI+reference mode, and SHALL NOT persist prompts, API keys, provider responses, document content, SMILES, or reasoning.

#### Scenario: Offline default and preference round-trip
- **GIVEN** an installation without a stored resolution preference
- **WHEN** Molecular Assistant loads its settings
- **THEN** the method defaults to AI+reference and external reference access defaults off
- **AND** an explicit external opt-in can be saved and restored without any secret or request data.
