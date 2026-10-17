# Molecular Depiction Visual Lab

## ADDED Requirements

### Requirement: Named experimental editorial presets
The lab SHALL provide a reusable named preset `ChemUSON Editorial AB3` with the exact owner-approved values for A1, A2, A3, B1, B2, and B3; A4 and B4 are not members of this calibration cohort. The preset SHALL only affect Working when explicitly selected and SHALL NOT change `CHEMDOODLE_LIKE`, production `DrawingStyle` defaults, Current, or the normal editor. The lab and command-line evidence export SHALL be able to load the preset without hand-editing JSON.

#### Scenario: Explicitly load editorial candidate
- **GIVEN** the lab is showing Current and Working with default parameters
- **WHEN** the user selects `ChemUSON Editorial AB3`
- **THEN** Working and its exports contain the exact approved parameter vector
- **AND** Current and global production defaults remain unchanged.

### Requirement: Editorial calibration eligibility is explicit
Every gallery case SHALL carry explicit editorial-calibration eligibility. Only A1, A2, A3, B1, B2, and B3 MAY inform the candidate. A4 and B4 SHALL be excluded from the approved calibration cohort. B4 SHALL remain renderable as a geometric stress fixture but SHALL be identified as non-editorial and excluded from preset selection and calibration summaries. Projection/diagnostic cases outside that six-case cohort SHALL NOT silently enter the AB3 selection corpus.

#### Scenario: Render B4 without editorial contamination
- **GIVEN** B4 remains in the gallery
- **WHEN** it is rendered or listed in an aggregate report/contact sheet
- **THEN** it is marked `geometric fixture`, `non-editorial`, and `excluded from preset selection / style calibration`
- **AND** editorial eligibility summaries omit it.
