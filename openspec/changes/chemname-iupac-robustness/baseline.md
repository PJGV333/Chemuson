# Baseline — ChemName / IUPAC robustness

Captured before production or test-data edits on **2026-10-08**.

## Repository identity

- Working directory: `/home/ccachyavgp/Documentos/ChemUSON`.
- Current branch: `chemname/iupac-robustness`, created from verified SHA `6aeef19028ffd047f8a39bc8f6063ea0b57210bf` (`release/v0.3.0-beta.1-prep`).
- Working tree was clean at capture (`git status --short`: empty).
- `release/v0.3.0-beta.1-prep` and `origin/release/v0.3.0-beta.1-prep` both resolve to `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`.
- `main` and `origin/main` both resolve to `f0fde603371255902bf0e630ca1ce032b8f16ad3`; no change is made to either ref.

## Reproduced source outputs

Using `smiles_to_molgraph` and `iupac_name` with `NameOptions(rdkit_isolated=False)` under an external 60-second timeout:

```text
CC(=O)N => 1-aminoethanamide
CC(=O)Cc1ccc(N)cc1 => 1-phenylpropan-2-one
CC(=O)Cc1cc(N)ccc1C => 1-phenylpropan-2-one
```

The last input loses both ring substituents. The analogous direct aryl ketone `CC(=O)c1cc(N)ccc1C` fails the current strict naming path with `MolecularViewNotSupported: Non-carbon in alkyl branch`; the safe default returns `N/D`. Direct aryl ketones remain out of scope.

The existing corpus contains 74 cases and expected the incorrect `^1-aminoethanamide$` for `CC(=O)N`; therefore the prior green corpus was not an independent correctness oracle for that case.

## Bounded baseline commands

- `timeout --signal=TERM --kill-after=5s 120s env PYTHONPATH=src pytest -q tests/test_chemname_*.py tests/test_iupac_ui.py`: **216 passed, 7 skipped in 33.07s**.
- `timeout --signal=TERM --kill-after=3s 60s python -m compileall -q src tests tools packaging`: **PASS**.
- `timeout --signal=TERM --kill-after=3s 30s ruff check src/chemuson/chemname tests/test_chemname_*.py --select F401,F811,F821,E722,E741`: **All checks passed**.
- A repository-wide `pytest --collect-only -q` had been started earlier in the campaign but its displayed output was truncated; no global collection count is asserted here.
- No monolithic `pytest -q` was run, in accordance with the inherited bounded-test/no-monolith constraint. No timed-out test was retried.

No production, test, corpus, roadmap, architecture catalog, package, release, or public-channel file had been changed at this baseline point.

## Stage 2/3 expansion baseline — 2026-10-08

Starting from the requested `a89e1f8562595db718a9d781b3183a4c2f5942de` on `chemname/iupac-robustness`; the working tree was clean and the branch matched `origin/chemname/iupac-robustness`. Release-prep and `main` refs matched their preserved SHAs above.

- Bounded existing ChemName/UI shard: `timeout --signal=TERM --kill-after=5s 120s env PYTHONPATH=src pytest -q tests/test_chemname_*.py tests/test_iupac_ui.py`: **226 passed, 7 skipped in 34.84s**.
- `timeout --signal=TERM --kill-after=3s 60s python -m compileall -q src tests tools packaging`: **PASS**.
- `timeout --signal=TERM --kill-after=3s 30s ruff check src/chemuson/chemname tests/test_chemname_*.py --select F401,F811,F821,E722,E741`: **PASS**.
- A structured PUG REST survey resolved **88 SMILES queries** to **78 unique PubChem CIDs**. For every unique structure, the returned PubChem connectivity SMILES canonicalized to the same RDKit graph as the input; formulas, registry names, and connectivity strings are recorded in `tests/data/chemname_iupac_reference_campaign.psv`.
- Stage 1 had **2** independently reference-backed structures. One (acetamide) overlaps this survey; the survey adds **77** unique structures. A follow-up PubChem-backed ester regression case (CID 12052398) is also recorded, bringing campaign-wide independent structure references to **80**.
- The original 78-row survey baseline classified **53 correct**, **16 incorrect**, **8 unsupported (`N/D`)**, and **1 reference pending**. The additive ester row was also incorrect at baseline, so the expanded 79-row fixture contains **53 correct**, **17 incorrect**, **8 unsupported**, and **1 reference pending**. `correct` includes 30 exact PubChem IUPACName matches and 23 documented valid systematic/retained variants; it does not mean every output is the PIN.
- Baseline family counts and overlap are in `corpus.md`. No production code was edited before these names and classifications were captured.
