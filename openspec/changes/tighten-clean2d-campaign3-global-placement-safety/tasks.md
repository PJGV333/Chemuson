# Tasks: Tighten Campaign 3 global placement safety

## Baseline and RED tests

- [x] Capture the current branch baseline and classify the four allowed failures plus the historical F401.
- [x] Add RED tests for explicit displacement-budget rejection, acceptance within budget, complete hard-gate evaluation, and deterministic safe competition.

## Safety implementation

- [x] Add a finite target/topology-dependent global displacement budget.
- [x] Replace the implicit displacement bypass with explicit global checks and JSON-safe evidence.
- [x] Keep all other hard gates strict and preserve the existing fallback.
- [x] Compare safe global/scaffold/unwrap candidates deterministically in complex preserve.

## Evidence and promotion

- [x] Update only the new corrective evidence with safety and competition fields.
- [x] Run the required focused tests, architecture, full suite, compileall, Ruff, and diff checks.
- [x] Validate this OpenSpec strictly and commit `Tighten Clean2D Campaign 3 global placement safety`.
- [ ] Archive the OpenSpec, repair canonical Purpose if needed, validate it, and commit `Archive Clean2D Campaign 3 safety closure`.
- [ ] Do not start Campaign 4 until this safety closure is promoted and archived.
