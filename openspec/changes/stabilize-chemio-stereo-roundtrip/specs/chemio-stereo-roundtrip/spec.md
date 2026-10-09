# Spec Delta — ChemIO stereo round-trip

## Purpose

Define preservation of molecular identity and explicitly specified stereochemistry across supported ChemIO SMILES and MOL/SDF conversions.

## ADDED Requirements

### Requirement: ChemIO conversions preserve molecular identity
ChemIO SHALL preserve molecular connectivity, atom and bond identity, bond orders, formal charges, isotopes, and supported atom/bond annotations across conversions. Serialization SHALL be deterministic under equivalent conditions. Tests SHALL compare chemical equivalence, not literal SMILES spelling.

#### Scenario: Round-trip a molecular graph
- **WHEN** a supported graph is serialized and re-imported through SMILES or MOL/SDF
- **THEN** atom/bond connectivity, element/isotope/charge, bond order, and molecular formula remain equivalent
- **AND** the serialization is deterministic under equal inputs/options.

### Requirement: Specified tetrahedral stereochemistry is preserved
A specified tetrahedral center SHALL remain chemically equivalent through each supported conversion, independent of SMILES neighbor order. Conversion MUST NOT turn specified stereo into unspecified/opposite stereo. A `@`/`@@` character comparison alone is insufficient.

#### Scenario: Opposite enantiomers round-trip
- **WHEN** either enantiomer of a supported chiral molecule is imported and exported
- **THEN** independent stereochemical descriptors identify the same enantiomer
- **AND** the opposite enantiomer remains chemically distinct.

### Requirement: Unspecified stereochemistry remains unspecified
ChemIO SHALL NOT assign tetrahedral or E/Z stereochemistry that is absent from the input. Visual wedge/hash presentation MUST NOT be invented to satisfy a display test.

#### Scenario: Potential center has no assignment
- **WHEN** an input contains a potentially stereogenic but unspecified center
- **THEN** the graph/output remains unspecified at that center and no artificial wedge/hash is added.

### Requirement: Supported E/Z annotations are retained
Explicit E/Z stereochemistry SHALL be preserved when the source and destination formats/backend support it. If specified E/Z cannot be represented faithfully, ChemIO SHALL fail explicitly or report an explicit warning rather than return a silently stereo-degraded structure.

#### Scenario: E/Z round-trip
- **WHEN** a supported E or Z alkene is exported and re-imported
- **THEN** RDKit independently assigns the same E/Z descriptor
- **AND** an unsupported conversion reports its limitation instead of erasing the descriptor.

### Requirement: MOL/SDF wedge parity is oriented
MOL/SDF stereo bond directions SHALL remain attached to the correct bond endpoint through internal parsing and serialization. Reordering atom/bond neighbors SHALL not reverse or erase tetrahedral identity.

#### Scenario: Reordered CTAB bond endpoints
- **WHEN** a MOL stereo bond is parsed with either endpoint first
- **THEN** its orientation/parity is retained and the effective stereoisomer is unchanged on export.

## Invariants

- Connectivity; atom/bond identity; bond orders; formal charges; isotopes.
- Specified tetrahedral centers and supported specified E/Z.
- No artificial assignment of unspecified stereo.
- Deterministic serialization for equivalent inputs/options.
- No change to Clean2D geometry or ChemName rules.
