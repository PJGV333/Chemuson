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

The machine-readable source of the original exact regression cases is `tests/data/chemname_acceptance_cases.yml`; each named case is exercised by `tools/chemname_acceptance.py`. The acceptance harness records status, exact name, input-build time, name-generation time, and total duration. Timing values are diagnostics, not hardware-sensitive gates.

## Stage 2/3 survey — 78 unique independently resolved structures

On **2026-10-08**, 88 PubChem PUG REST queries were deduplicated to 78 unique molecular structures. Each row in [`tests/data/chemname_iupac_reference_campaign.psv`](../../../tests/data/chemname_iupac_reference_campaign.psv) records the submitted SMILES, formula, PubChem CID, PubChem `IUPACName`, returned connectivity SMILES, exact pre-fix ChemName output, and baseline classification. All 78 input/connectivity pairs produced the same RDKit-canonicalized graph. The 88 requests included overlapping representations/cross-family cases; metrics below count unique structures within a family and are intentionally not additive.

Stage 1 already had two independently referenced structures. Acetamide is one of these 78; the substituted aryl ketone is not in this survey. Thus the survey adds **77 new structures**, for **79 unique independently referenced structures** across the campaign. PubChem `IUPACName` is an independent registry result, not assumed automatically to be a PIN. The baseline assessed 53 outputs as chemically complete/correct, including **30 exact PubChem name matches and 23 valid systematic/retained-name variants**; 16 as incorrect/incomplete or contrary to the applicable preferred suffix convention; 8 as unsupported (`N/D`); and 1 as reference-pending.

### Baseline classifications by overlapping family

| Family | Structures | Correct | Incorrect | Unsupported | Reference pending |
|---|---:|---:|---:|---:|---:|
| Carboxylic acids (including multifunctional acids) | 12 | 8 | 2 | 2 | 0 |
| Esters | 9 | 0 | 5 | 3 | 1 |
| Aldehydes | 8 | 6 | 1 | 1 | 0 |
| Ketones | 8 | 6 | 1 | 1 | 0 |
| Alcohols | 8 | 6 | 2 | 0 | 0 |
| Amines | 8 | 7 | 1 | 0 | 0 |
| Amides | 8 | 5 | 3 | 0 | 0 |
| Aromatic structures | 15 | 14 | 0 | 1 | 0 |
| Multifunctional structures | 12 | 10 | 2 | 0 | 0 |

A family tag can overlap another tag (for example, an amino acid is both acid and multifunctional). `correct` accepts a documented valid systematic or retained variant even if it is not PubChem's exact string or the PIN. The single pending adjudication is whether `1-acetoxyethane` is acceptable general prefix-mode nomenclature for ethyl acetate; the reference-backed preferred form is `ethyl acetate` (Blue Book P-65.6.3.3.1).

### Reproducible normative anchors

- Seniority and selection of suffix classes: IUPAC Blue Book (2013), P-41 and P-44 ([P-4](https://iupac.qmul.ac.uk/BlueBook/P4.html), [P-5](https://iupac.qmul.ac.uk/BlueBook/P5.html)).
- Amines: P-62; hydroxy compounds and multiplicative `diol`: P-63.1.2; ketones, `oxo`, and polyfunctional ketones: P-64, especially P-64.2.1.2/P-64.7.
- Carboxylic acids and polyacids: P-65.1, especially P-65.1.2.2; ester organyl-plus-anion names: P-65.6.3.3.1; see [P-6](https://iupac.qmul.ac.uk/BlueBook/P6.html).
- N-substituted and multiple amides: P-66.1.1.1.1.1/P-66.1.1.3.1; mono-/dialdehydes: P-66.6.1.1.1; see [P-6a](https://iupac.qmul.ac.uk/BlueBook/P6a.html).

### Defects and unsupported candidates captured before fixes

- Acid parent selection/suffix: `CC(C)C(=O)O` and benzoic acid return `N/D`; `O=C(O)C(=O)O` becomes `2-carboxyethanoic acid` (adds a carbon), and succinic acid `O=C(O)CCC(=O)O` becomes `4-carboxybutanoic acid` (also adds a carboxyl carbon).
- Ester component loss: methyl/ethyl propanoate, methyl acetate, methyl propenoate, and ethyl 2-hydroxypropanoate return only the acid-derived `...oate` name; aromatic benzoates and methyl 2-methylpropanoate return `N/D`.
- Branched aldehyde `CC(C)C=O` returns `N/D`; propanedial returns `3-oxopropanal` instead of the documented preferred dial suffix.
- Same-class suffixes are not combined for `OCCO`, `OCCCO`, `NCCN`, and `CC(=O)CCC(=O)C`; names instead demote one identical function to `hydroxy`/`amino`/`oxo`.
- `CNC(C)=O`, `CCNC(C)=O`, and `NC(=O)CC(=O)N` omit or misrepresent amide N-substitution/multiplicity. `CCOC(=O)CCN` omits the ethyl ester component.
- Acetophenone and 4-hydroxyacetophenone return `N/D`; their simple structures have independent PubChem references and are stage-3 support candidates, not deliberate boundaries.
