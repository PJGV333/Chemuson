## Context

Clean 2D baseline diff review summaries are test artifacts that describe observed differences in baseline reports. They currently help reviewers see that something changed, but they do not provide a stable severity recommendation for how a human should treat each difference.

This design keeps the policy in test-only infrastructure. It consumes existing review summaries/items and produces a structured advisory decision for reviewers. The decision is intentionally observational and must not fail CI or change production Clean 2D behavior.

## Goals / Non-Goals

**Goals:**
- Provide a small, deterministic helper for classifying Clean 2D diff review summaries.
- Preserve visibility into every reason, including known-delicate cases.
- Make geometry-risk review-required while marking it observational-only.
- Keep output JSON-stable for tests and reviewer tooling.

**Non-Goals:**
- Add automatic gates, CI failures, or enforcement.
- Change Clean 2D layout, ranking, selection, UI, backends, public document format, RDKit/CoordGen behavior, or production code.
- Repair, improve, or normalize geometry.
- Refactor Clean 2D engine code.

## Decisions

- Implement the policy in `tests/clean2d_regression/blocking_policy.py` rather than production code. This keeps the policy available to tests and reviewers without creating a runtime dependency or enforcement surface.
- Expose one evaluator that accepts an existing review summary mapping. This avoids changing the baseline diff review summary format and keeps the new policy layered on top of current review artifacts.
- Return a frozen dataclass with primitive JSON-friendly fields plus a `to_dict()` method. This gives tests a stable object model while preserving deterministic serialization with `json.dumps(..., allow_nan=False, sort_keys=True)`.
- Treat `known-delicate-change` as visible and allowed only when it is the highest severity present. Mixed summaries still surface stronger risks, so known-delicate evidence cannot hide contract-risk items.
- Treat geometry-risk as `review-required` and `observational_only=True`. The policy recommends review but explicitly remains non-gating.

## Risks / Trade-offs

- Policy labels may be mistaken for enforcement. Mitigation: include explicit `observational_only` output and tests that verify no gate behavior is introduced.
- Existing review item shapes may vary. Mitigation: accept mapping-like summaries and inspect stable fields such as `decision`, `category`, `severity`, `reason`, and `known_delicate` without mutating input.
- Future gates may need stronger metadata. Mitigation: keep `future-gate-candidate` documented as a severity vocabulary item but do not enforce it in this change.
