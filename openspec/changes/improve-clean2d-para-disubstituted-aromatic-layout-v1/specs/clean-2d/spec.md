## ADDED Requirements

### Requirement: Para-disubstituted aromatic template
Clean2D SHALL provide a deterministic safe template candidate for a complete, neutral, stereo-free carbocyclic aromatic six-ring with exactly two opposite allowed one-atom terminal substituents.

#### Scenario: Para-disubstituted benzene
- **WHEN** Quick or Publication Clean2D evaluates a matching para-disubstituted benzene
- **THEN** it MAY select `para_disubstituted_aromatic_template` with both substituents radial and opposite.

#### Scenario: Unsupported substitution pattern
- **WHEN** the structure is ortho, meta, fused, charged, stereo-marked, partially selected, or has a non-terminal substituent
- **THEN** Clean2D SHALL NOT generate the para-disubstituted aromatic template.

### Requirement: Clean2D visual score direction
`visual_score` SHALL be interpreted as lower-is-better when comparing Clean2D candidate quality.

#### Scenario: Compare candidate quality
- **WHEN** two candidates have different visual scores
- **THEN** the candidate with the lower visual score has the better visual score
