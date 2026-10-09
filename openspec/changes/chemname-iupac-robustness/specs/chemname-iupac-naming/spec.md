## MODIFIED Requirements

### Requirement: Primary amides are named as amides

ChemName SHALL recognize a primary carboxamide as the principal amide function and SHALL NOT also name its nitrogen as an independent amino/amine substituent. For the exact structure `CC(=O)N`, the requested systematic output SHALL be `ethanamide` (the Blue Book retained PIN `acetamide` remains an accepted nomenclature variant generally, but is not this campaign's display expectation).

#### Scenario: Acetamide systematic output
- **GIVEN** the exact neutral structure `CC(=O)N`
- **WHEN** ChemName names it using the default supported path
- **THEN** the output is exactly `ethanamide`
- **AND** it contains no `amino` prefix or duplicated N-derived group

### Requirement: Supported aryl substituents retain their decorations

When ChemName names a neutral, singly attached benzene as a phenyl substituent, it SHALL preserve every recognized supported substituent and its attachment-relative locant. A substituted phenyl name SHALL be grouped unambiguously in the enclosing name. If the ring connection, substituent metadata, or stereochemistry cannot be represented by this naming path, ChemName SHALL fail closed rather than return an incomplete bare `phenyl` name.

#### Scenario: Multisubstituted phenyl on a ketone chain
- **GIVEN** the exact neutral structure `CC(=O)Cc1cc(N)ccc1C`
- **WHEN** ChemName names the propan-2-one parent with the substituted ring
- **THEN** the output is exactly `1-(5-amino-2-methylphenyl)propan-2-one`
- **AND** both the amino and methyl substituents and their locants are retained

#### Scenario: Unsupported stereochemistry in an aryl substituent
- **GIVEN** a phenyl substituent carrying stereochemical metadata on a branch that this naming path cannot describe
- **WHEN** ChemName attempts to name the structure
- **THEN** it returns `N/D` under the default safe option (or raises `ChemNameNotSupported` when `return_nd_on_fail` is disabled)
- **AND** it SHALL NOT discard the stereochemistry and return a bare or partially decorated phenyl name

#### Scenario: Unsupported aryl connectivity or decoration
- **GIVEN** a ring with multiple connections to the selected parent, a cross-link, a charged/isotopically modified unsupported decoration, or an unrecognized substituent
- **WHEN** ChemName attempts to use it as a simple phenyl substituent
- **THEN** it fails closed rather than emit a name that omits connectivity, functional groups, isotope/charge state, or stereo
