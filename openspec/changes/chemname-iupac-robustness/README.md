# ChemName / IUPAC robustness

Bounded correctness campaign created on branch `chemname/iupac-robustness` from verified source SHA `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`.

The first implementation tranche corrected two reproducible behaviors:

1. `CC(=O)N` is an amide and must not be decorated as `1-aminoethanamide`. The requested systematic display is `ethanamide`; the OpenSpec records that the IUPAC Blue Book also recognizes retained `acetamide` as the PIN, so the campaign does not mislabel the output convention.
2. `CC(=O)Cc1cc(N)ccc1C` must not silently lose the amino and methyl substituents when its ring is named as a phenyl substituent. The independently retrieved PubChem IUPAC name is `1-(5-amino-2-methylphenyl)propan-2-one` (CID 118802021).

An expanded campaign was authorized on 2026-10-08 without resetting this OpenSpec or its prior commits. It adds a measured 78-structure survey plus a separately verified ester regression (79 fixture rows; 80 unique campaign references), with staged work on functional-group priority, aromatic parents, structure accounting, and determinism; see the appended tasks and corpus/validation evidence. This remains bounded ChemName work, not a claim of general IUPAC compliance. Formula/isotope/charge presentation, RDKit isolation, GUI behavior, Clean2D, persistence, packaging, release metadata, and public channels remain outside scope.
