# Reference-backed ChemName regression corpus

The original survey records were retrieved on **2026-10-08**; the methylammonium reference was retrieved on **2026-10-09**. Names are compared against exact input structures; no ChemName-generated value is used as its own reference.

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

## C3 — Methylammonium charge preservation

| Field | Value |
|---|---|
| Isomeric SMILES | `C[NH3+]` |
| Formula | `CH6N+` |
| Baseline ChemName | `methan-1-amine` (loses the cation) |
| Required campaign output | `methylazanium` |
| Independent registry | PubChem CID 644041; `IUPACName`: `methylazanium`; connectivity `C[NH3+]` |
| Retrieval date | 2026-10-09 |

The exact cation now has a dedicated positive assertion. The adjacent ethylammonium form remains outside this narrowly referenced support case and fails closed rather than inheriting the methyl name.

Reference:
- [PubChem CID 644041](https://pubchem.ncbi.nlm.nih.gov/compound/644041)
- [PubChem PUG REST properties for `C[NH3+]`](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/C%5BNH3%2B%5D/property/MolecularFormula,IUPACName,ConnectivitySMILES/JSON)

## Negative / safety sentinels

- `CC(=O)Cc1cc(N)ccc1[C@H](C)O` carries a stereogenic hydroxyethyl branch on the phenyl group. The isolated SMILES import represents the stereochemistry as hashed-bond metadata; the campaign's exact safe expectation is `N/D` until this path can emit and reference the descriptor. A separate graph-level unit test sets `stereo_cip=R` to cover atom metadata as well.
- Unsupported ring-to-parent multiple connections, cross-links, or decorations whose charge/isotope/function/stereo is not represented must likewise fail closed. Unit sentinels exercise charged/isotopic branch metadata and multiple ring-to-parent attachments.
- `CC(=O)c1cc(N)ccc1C` currently fails the engine's supported path and returns `N/D` by default. Direct aryl ketones remain outside this tranche; the change must not turn this unsupported case into a partial name.

## Reproducibility and metric policy

The machine-readable source of the original exact regression cases is `tests/data/chemname_acceptance_cases.yml`; each named case is exercised by `tools/chemname_acceptance.py`. The acceptance harness records status, exact name, input-build time, name-generation time, and total duration. Timing values are diagnostics, not hardware-sensitive gates.

## Stage 2/3 survey — 78 unique independently resolved structures

On **2026-10-08**, 88 PubChem PUG REST queries were deduplicated to 78 unique molecular structures. The 80-row [`tests/data/chemname_iupac_reference_campaign.psv`](../../../tests/data/chemname_iupac_reference_campaign.psv) contains those 78 survey structures, one independently verified follow-up regression case (CID 12052398), and the separately resolved methylammonium charge case (CID 644041). Each row records the submitted SMILES, formula, PubChem CID, PubChem `IUPACName`, returned connectivity SMILES, exact pre-fix ChemName output, baseline classification, and final target. All input/connectivity pairs produced the same RDKit-canonicalized graph. The 88 survey requests included overlapping representations/cross-family cases; metrics below count unique structures within a family and are intentionally not additive.

Stage 1 already had two independently referenced structures. Acetamide is one of the 78 survey structures; the substituted aryl ketone is separate. Thus the survey adds **77 new structures**, for **79** at the initial campaign checkpoint; the follow-up ester brings that count to **80**, and methylammonium to **81 total independently referenced campaign structures**. PubChem `IUPACName` is an independent registry result, not assumed automatically to be a PIN. The original 78-row baseline assessed 53 outputs as chemically complete/correct, including **30 exact PubChem name matches and 23 valid systematic/retained-name variants**; 16 as incorrect/incomplete; 8 as unsupported (`N/D`); and 1 as reference-pending. The follow-up ester and methylammonium cases were also incorrect at baseline, so the 80-row fixture baseline is **53 correct, 18 incorrect, 8 unsupported, and 1 reference-pending**.

#### Follow-up ester reference

The existing amino-oxo ester regression graph in `tests/test_chemname_pr27.py` corresponds to `COC(=O)C(N)C(=O)C`, formula `C5H9NO3`. PubChem PUG REST resolves CID 12052398, connectivity `CC(=O)C(C(=O)OC)N`, and `IUPACName` **methyl 2-amino-3-oxobutanoate**. ChemName's previous `2-amino-3-oxobutanoate` omitted the methyl ester component. This independent identity check backs the corrected complete name and also served as a negative control for the earlier family regression.

### Final measured fixture coverage after family fixes

A strict final audit of all 80 fixture structures (`rdkit_isolated=False`, `return_nd_on_fail=False`) produced no exceptions or `N/D`: 78 outputs match their recorded target strings and the remaining two are the documented systematic N-substituted amide forms `N-methylethanamide` and `N-ethylethanamide`. The formerly pending ethyl acetate case now matches its reference target; methylammonium is named `methylazanium`, matching PubChem CID 644041 and preserving the cation. The separate stage-1 substituted-aryl reference also remains covered by its exact regression test, bringing campaign-wide independent references to 81. These measurements describe only this finite corpus; they do not establish general IUPAC compliance. The per-structure final audit is enforced by `tests/test_chemname_reference_campaign.py`.

### Baseline classifications by overlapping family

| Family | Structures | Correct | Incorrect | Unsupported | Reference pending |
|---|---:|---:|---:|---:|---:|
| Carboxylic acids (including multifunctional acids) | 12 | 8 | 2 | 2 | 0 |
| Esters | 10 | 0 | 6 | 3 | 1 |
| Aldehydes | 8 | 6 | 1 | 1 | 0 |
| Ketones | 8 | 6 | 1 | 1 | 0 |
| Alcohols | 8 | 6 | 2 | 0 | 0 |
| Amines | 9 | 7 | 2 | 0 | 0 |
| Amides | 8 | 5 | 3 | 0 | 0 |
| Aromatic structures | 15 | 14 | 0 | 1 | 0 |
| Multifunctional structures | 13 | 10 | 3 | 0 | 0 |

A family tag can overlap another tag (for example, an amino acid is both acid and multifunctional). `correct` accepts a documented valid systematic or retained variant even if it is not PubChem's exact string or the PIN. The one `reference_pending` classification is historical: the baseline output `1-acetoxyethane` had not been adjudicated. The current output is the independently backed `ethyl acetate` (CID 8857; Blue Book P-65.6.3.3.1), so no current corpus case remains pending.

#### Current final status by overlapping family

| Family | Structures | Exact target strings | Adjudicated systematic variants | Unsupported/other |
|---|---:|---:|---:|---:|
| Carboxylic acids | 12 | 12 | 0 | 0 |
| Esters | 10 | 10 | 0 | 0 |
| Aldehydes | 8 | 8 | 0 | 0 |
| Ketones | 8 | 8 | 0 | 0 |
| Alcohols | 8 | 8 | 0 | 0 |
| Amines | 9 | 9 | 0 | 0 |
| Amides | 8 | 6 | 2 | 0 |
| Aromatic structures | 15 | 15 | 0 | 0 |
| Multifunctional structures | 13 | 13 | 0 | 0 |

Family rows overlap by design and therefore do not sum to the 80 unique fixture structures. The two systematic variants are the N-methyl and N-ethyl ethanamide outputs; PubChem's retained acetamide strings remain recorded as independent registry names.

### Reproducible normative anchors

- Seniority and selection of suffix classes: IUPAC Blue Book (2013), P-41 and P-44 ([P-4](https://iupac.qmul.ac.uk/BlueBook/P4.html), [P-5](https://iupac.qmul.ac.uk/BlueBook/P5.html)).
- Amines: P-62; hydroxy compounds and multiplicative `diol`: P-63.1.2; ketones, `oxo`, and polyfunctional ketones: P-64, especially P-64.2.1.2/P-64.7.
- Carboxylic acids and polyacids: P-65.1, especially P-65.1.2.2; ester organyl-plus-anion names: P-65.6.3.3.1; see [P-6](https://iupac.qmul.ac.uk/BlueBook/P6.html).
- N-substituted and multiple amides: P-66.1.1.1.1.1/P-66.1.1.3.1; mono-/dialdehydes: P-66.6.1.1.1; see [P-6a](https://iupac.qmul.ac.uk/BlueBook/P6a.html).

### Name-form adjudications used during fixes

For `CNC(C)=O` and `CCNC(C)=O`, PubChem records `N-methylacetamide` and `N-ethylacetamide`; the implementation now emits the fully specified systematic-parent forms `N-methylethanamide` and `N-ethylethanamide`. The corpus retains PubChem strings as independent registry references, not automatic exact-output or PIN oracles. This follows the campaign's recorded distinction between systematic `ethanamide` and retained `acetamide`; both forms preserve the N-substituent and parent connectivity.

### Stage 3 aromatic cases now supported

The exact structures already in the corpus now produce their recorded target names: benzoic acid (CID 243), methyl benzoate (CID 7150), ethyl benzoate (CID 7165), acetophenone (CID 7410), and 4-hydroxyacetophenone (CID 7469). Acetophenone outputs use explicit locant `1` (`1-phenylethan-1-one`), whereas PubChem's registry string is `1-phenylethanone`; the corpus target records the chosen systematic form. Fixed atom-order pairs cover both aromatic carboxyl parents and aryl ketones. Direct aryl-ketone support is intentionally allowlisted to acetophenone and its para-hydroxy case; other substituted direct aryl ketones such as `CC(=O)c1cc(N)ccc1C` remain `N/D` even though the separate phenyl-substituent path can name `CC(=O)Cc1cc(N)ccc1C`. Fused, charged, isotopic, and stereogenic direct aryl-ketone decorations remain `N/D`; phenylacetic acid remains distinct from benzoic acid.

### Defects and unsupported candidates captured before fixes

- Acid parent selection/suffix: `CC(C)C(=O)O` and benzoic acid return `N/D`; `O=C(O)C(=O)O` becomes `2-carboxyethanoic acid` (adds a carbon), and succinic acid `O=C(O)CCC(=O)O` becomes `4-carboxybutanoic acid` (also adds a carboxyl carbon).
- Ester component loss: methyl/ethyl propanoate, methyl acetate, methyl propenoate, ethyl 2-hydroxypropanoate, ethyl 3-aminopropanoate, and methyl 2-amino-3-oxobutanoate return only the acid-derived `...oate` name; aromatic benzoates and methyl 2-methylpropanoate return `N/D`.
- Branched aldehyde `CC(C)C=O` returns `N/D`; propanedial returns `3-oxopropanal` instead of the documented preferred dial suffix.
- Same-class suffixes are not combined for `OCCO`, `OCCCO`, `NCCN`, and `CC(=O)CCC(=O)C`; names instead demote one identical function to `hydroxy`/`amino`/`oxo`.
- `CNC(C)=O`, `CCNC(C)=O`, and `NC(=O)CC(=O)N` omit or misrepresent amide N-substitution/multiplicity. `CCOC(=O)CCN` omits the ethyl ester component. `C[NH3+]` was named as neutral `methan-1-amine`, losing charge; its reference-backed final name is `methylazanium` (PubChem CID 644041).
- Acetophenone and 4-hydroxyacetophenone return `N/D`; their simple structures have independent PubChem references and are stage-3 support candidates, not deliberate boundaries.
