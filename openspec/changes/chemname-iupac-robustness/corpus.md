# Reference-backed ChemName regression corpus

Retrieval date for all online records: **2026-10-08**. Names are compared against exact input structures; no ChemName-generated value is used as its own reference.

## C1 — Primary carboxamide

| Field | Value |
|---|---|
| ID | `smiles_acetamide` (existing acceptance case) |
| Isomeric SMILES | `CC(=O)N` |
| Baseline ChemName | `1-aminoethanamide` |
| Required campaign output | `ethanamide` |
| Formula/connectivity | CH3–C(=O)–NH2; one carboxamide, not a separate amine plus amide |
| Normative rule | IUPAC Blue Book (2013), P-66.1.1.1.1.1: substitutive carboxamide names are formed by adding `amide` to the parent hydride, eliding final `e` before `a`; the parent hydride is ethane. P-66.1.1.1.2.1 independently lists retained `acetamide` as the PIN. |
| Independent registry cross-check | PubChem CID 178, PUG REST `IUPACName`: `acetamide` (the retained PIN; confirms structure and nomenclature variant, not the campaign's requested systematic display). |

**Interpretation:** The exact requested output `ethanamide` is the systematic substitutive form under the cited rule. It is not described here as the Blue Book PIN; that PIN is `acetamide`. The actual defect is unambiguous: `1-aminoethanamide` wrongly treats the carboxamide nitrogen as an independent amino substituent.

References:
- [IUPAC Blue Book, P-66 (amides)](https://iupac.qmul.ac.uk/BlueBook/P6a.html#66)
- [IUPAC Blue Book, P-2 (parent hydrides / ethane)](https://iupac.qmul.ac.uk/BlueBook/P2.html)
- [PubChem CID 178](https://pubchem.ncbi.nlm.nih.gov/compound/178)
- [PubChem PUG REST query for `CC(=O)N`](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CC%28%3DO%29N/property/IUPACName/JSON)

## C2 — Multisubstituted phenyl on a ketone chain

| Field | Value |
|---|---|
| ID | `smiles_substituted_aryl_ketone` |
| Isomeric SMILES | `CC(=O)Cc1cc(N)ccc1C` |
| Baseline ChemName | `1-phenylpropan-2-one` (drops both amino and methyl substituents) |
| Required campaign output | `1-(5-amino-2-methylphenyl)propan-2-one` |
| Independent registry | PubChem CID 118802021, `IUPACName`: `1-(5-amino-2-methylphenyl)propan-2-one`; PUG REST also returned connectivity SMILES `CC1=C(C=C(C=C1)N)CC(=O)C`. |
| Scope | A neutral, single-attachment phenyl substituent with supported amino and methyl decorations. This does not claim support for direct acyl-substituted benzene (e.g. acetophenone), fused rings, nested aryl expansion, or stereochemical branches. |

Reference:
- [PubChem CID 118802021](https://pubchem.ncbi.nlm.nih.gov/compound/118802021)
- [PubChem PUG REST query for the exact input SMILES](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CC%28%3DO%29Cc1cc%28N%29ccc1C/property/IUPACName,CanonicalSMILES/JSON)

## Negative / safety sentinels

- `CC(=O)Cc1cc(N)ccc1[C@H](C)O` carries a stereogenic hydroxyethyl branch on the phenyl group. The isolated SMILES import represents the stereochemistry as hashed-bond metadata; the campaign's exact safe expectation is `N/D` until this path can emit and reference the descriptor. A separate graph-level unit test sets `stereo_cip=R` to cover atom metadata as well.
- Unsupported ring-to-parent multiple connections, cross-links, or decorations whose charge/isotope/function/stereo is not represented must likewise fail closed. Unit sentinels exercise charged/isotopic branch metadata and multiple ring-to-parent attachments.
- `CC(=O)c1cc(N)ccc1C` currently fails the engine's supported path and returns `N/D` by default. Direct aryl ketones remain outside this tranche; the change must not turn this unsupported case into a partial name.

## Reproducibility and metric policy

The machine-readable source of exact regression cases is `tests/data/chemname_acceptance_cases.yml`; each named case is exercised by `tools/chemname_acceptance.py`. The acceptance harness records status, exact name, input-build time, name-generation time, and total duration. The gating metric is **exact-name correctness / no unexpected fail or error**; timing values are reported diagnostically and no hardware-sensitive latency threshold is invented. Unsupported sentinels gate on explicit `N/D`/not-supported behavior.
