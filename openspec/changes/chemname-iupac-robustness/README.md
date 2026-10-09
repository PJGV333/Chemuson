# ChemName / IUPAC robustness

Bounded correctness campaign created on branch `chemname/iupac-robustness` from verified source SHA `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`.

This tranche addresses only two reproducible behaviors:

1. `CC(=O)N` is an amide and must not be decorated as `1-aminoethanamide`. The requested systematic display is `ethanamide`; the OpenSpec records that the IUPAC Blue Book also recognizes retained `acetamide` as the PIN, so the campaign does not mislabel the output convention.
2. `CC(=O)Cc1cc(N)ccc1C` must not silently lose the amino and methyl substituents when its ring is named as a phenyl substituent. The independently retrieved PubChem IUPAC name is `1-(5-amino-2-methylphenyl)propan-2-one` (CID 118802021).

See [`corpus.md`](corpus.md) for exact structures, citations, baseline outputs and scope. See [`baseline.md`](baseline.md), [`design.md`](design.md), [`tasks.md`](tasks.md), [`specs/chemname-iupac-naming/spec.md`](specs/chemname-iupac-naming/spec.md), and [`validation.md`](validation.md) for the contract and evidence. The campaign does not change other nomenclature rules, formula/isotope/charge presentation, RDKit isolation, GUI behavior, Clean2D, persistence, release metadata, or public channels.
