# canvas-selection-clipboard-policy Specification

## Purpose
TBD - created by archiving change extract-canvas-selection-clipboard-policy. Update Purpose after archive.
## Requirements
### Requirement: Clipboard policy and codec preserve formats

The selection clipboard module SHALL expose the existing MIME strings, pasteable
format policy, and UTF-8 JSON codec without owning clipboard, scene or undo
operations.

#### Scenario: supported formats
- **GIVEN** MIME data with any existing ChemUSON, Molfile, URL, text or image
  format
- **WHEN** pasteability is queried
- **THEN** the result is true

#### Scenario: payload roundtrip
- **GIVEN** a selection payload dictionary
- **WHEN** it is encoded and decoded
- **THEN** keys and values are preserved

### Requirement: Large-selection and bond-copy policies remain deterministic

The module SHALL preserve the existing atom/bond thresholds and duplicate-bond
priority behavior without mutating a model or scene.

#### Scenario: large selection threshold
- **GIVEN** counts at or above either configured threshold
- **WHEN** large-selection policy runs
- **THEN** it returns true

