# Validation — molecular identity verification

## Checks

- Identity verification unit tests: **6 passed**, fully offline, including equivalent, mismatch, missing reference, resolver error and non-applicable prompts.
- `cholesterol-semantic-mismatch-01` is included in evaluator data and verified offline; evaluator metadata records identity status and reference identifier without changing layout acceptance.
- Assistant controller/UI coverage: identity verification runs in the worker; separate UI label distinguishes ChemIO validity from requested identity; mismatch insertion requires a second explicit confirmation.
- `tests/architecture`: **278 passed** after cataloging the M10→M16 runtime edge.
- Strict OpenSpec: **52 passed, 0 failed**; changed-file selected Ruff: **PASS**; `git diff --check`: **PASS**.
- No direct RDKit import, network-based test resolver, or saved identity result was introduced. Canonicalization uses the existing isolated ChemIO worker.

The live Qwen smoke used an offline injected reference resolver only. Successful ethanol/caffeine outputs were verified against offline fixtures; model failures did not enter identity verification. Full suite limitations and unrelated chemistry failures are in `../stabilize-ai-molecular-assistant-integration/validation.md`.
