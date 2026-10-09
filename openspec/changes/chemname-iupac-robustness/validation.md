# Validation — ChemName / IUPAC robustness

## Baseline before implementation

- Branch/base and clean status: [`baseline.md`](baseline.md).
- Bounded ChemName/UI baseline: **216 passed, 7 skipped in 33.07s**.
- Compileall and scoped ChemName Ruff: **PASS**.
- Exact defective outputs and independent references: [`corpus.md`](corpus.md).

## Red reproduction

Before changing production logic:

- `timeout 60s env PYTHONPATH=src pytest -q tests/test_chemname_iupac_robustness.py`: **3 expected failures**. The observations were `1-aminoethanamide`, `1-phenylpropan-2-one`, and a stereo-bearing phenyl decoration incorrectly flattened to `1-phenylpropan-2-one`.
- The pre-fix acceptance subset reported **1 passed, 3 failed, 0 errors**: the two missing-substituent/amide expectations failed; the direct acyl-substituted benzene sentinel correctly remained `N/D`.

These were deliberate red tests, not final unresolved failures.

## Post-change correctness and metrics

- Exact bounded acceptance subset (`smiles_acetamide`, both atom-order forms of the substituted aryl ketone, direct aryl ketone unsupported, stereogenic aryl branch unsupported), under an external 60-second timeout: **5 passed, 0 failed, 0 skipped, 0 errors**. Exact outputs:
  - `CC(=O)N` → `ethanamide`.
  - `CC(=O)Cc1cc(N)ccc1C` → `1-(5-amino-2-methylphenyl)propan-2-one`.
  - `Cc1ccc(N)cc1CC(=O)C` → `1-(5-amino-2-methylphenyl)propan-2-one`.
  - `CC(=O)c1cc(N)ccc1C` → `N/D` (direct acyl-substituted benzene remains out of scope).
  - `CC(=O)Cc1cc(N)ccc1[C@H](C)O` → `N/D` (unsupported stereo is not discarded).
- Acceptance harness per-case total durations on this host: **379.19 ms**, **357.34 ms**, **353.68 ms**, **141.83 ms**, **142.76 ms**, respectively. These are diagnostic values, not machine-sensitive thresholds.
- Final bounded ChemName/UI shard, `timeout 120s env PYTHONPATH=src pytest -q tests/test_chemname_*.py tests/test_iupac_ui.py`: **226 passed, 7 skipped in 35.73s**. This includes the full 78-case acceptance dataset and packaged source smoke.
- Architecture suite, `timeout 60s env PYTHONPATH=src pytest -q tests/architecture`: **280 passed in 9.84s**.
- `timeout 60s python -m compileall -q src tests tools packaging`: **PASS**.
- Scoped Ruff (`F401,F811,F821,E722,E741`) over changed ChemName modules/tests: **PASS**.
- `openspec validate chemname-iupac-robustness --strict`: **PASS**.
- Acceptance JSON parse and `git diff --check`: **PASS**.

An additional unscoped default Ruff scan is non-clean with **58 legacy style diagnostics** across the touched modules (mostly pre-existing import ordering, old typing aliases, broad catches and modernization rules). The added phenyl helper has no remaining diagnostic in that scan. No unrelated style cleanup was made; the repository-mandated scoped rules pass.

## Boundaries and delivery

- No monolithic `pytest -q` was run. No test timed out; no timeout was retried.
- No Clean2D, ChemIO/RDKit worker, GUI, `.cmsn`, architecture catalog, dependency, release workflow, or public-channel change was made.
- `ethanamide` is the requested systematic display form; the Blue Book PIN is retained `acetamide`, as documented in `corpus.md`.
- Remaining unsupported cases fail closed; this tranche does not assert general IUPAC compliance.
- Local commits: `d7823bd` (scope/baseline) and `bdfd00233b6f1c7c8c3e3c67a65f83acb1dd4cde` (ChemName implementation/regressions). The branch is `chemname/iupac-robustness`; `release/v0.3.0-beta.1-prep` and its origin still resolve to `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`; `main` and `origin/main` remain `f0fde603371255902bf0e630ca1ce032b8f16ad3`. No push, merge, tag, release, or public-channel operation occurred.
