## 1. Schema & Helpers (Foundation)

- [x] 1.1 Create `clean_2d_quality` package with __init__.py and lightweight contract definitions using dataclass, TypedDict, plain dict helpers, or manual validation; do not add pydantic or another new dependency unless it already exists in the project.
- [x] 1.2 Implement helper function to construct diagnostics from raw data, enforcing field constraints (state enum, reason enum when required, score range 0–1, non-empty source).
- [x] 1.3 Add validation logic that raises clear errors if contract is violated (e.g., missing reason for rejected or failed-controlled state, out-of-range reporting score).

## 2. Contract Tests (Before Code Changes)

- [x] 2.1 Write unit tests validating the schema: correct field types, enum values, and required/optional semantics per spec requirements.
- [x] 2.2 Write contract test for successful diagnostic: GIVEN a molecule passes all checks, WHEN diagnostic is generated, THEN state=applied, no reason, score in (0,1], source non-empty.
- [x] 2.3 Write contract test for rejected diagnostic: GIVEN a molecule fails safety with stereo-risk, WHEN diagnostic is generated, THEN state=rejected, reason=stereo-risk, low score near 0, source=safety.
- [x] 2.4 Write contract test verifying internal metadata preservation without affecting stable contract fields.
- [x] 2.5 Write contract test validating score normalization for reporting only: GIVEN raw metrics from different modules, WHEN diagnostic is generated, THEN score is always in [0,1] regardless of original scale.
- [x] 2.6 Write contract tests for reason rules: rejected and failed-controlled require reason; no-op and preserve-only allow optional reason; applied omits reason unless an explicit non-contract diagnostic path includes it.

## 3. Minimal Module Adapters

- [x] 3.1 Identify existing clear diagnostic surfaces in engine.py, safety.py, depiction_candidates.py, and Clean2DController without refactoring their internals.
- [x] 3.2 Add minimal wrapper for engine.py only if an existing diagnostic surface can be adapted without changing layout behavior.
- [x] 3.3 Add minimal wrapper for safety.py only if current checks can be mapped to stable reasons without changing safety decisions.
- [x] 3.4 Add minimal wrapper for depiction_candidates.py only for reporting candidate source and reporting score; do not change ranking or candidate selection.
- [x] 3.5 Add minimal Clean2DController pass-through only if needed to expose stable diagnostics; do not change UI messages except minimally to connect diagnostics.
- [x] 3.6 Record any non-trivial duplication or broad migration need as follow-up instead of removing it in this change.

## 4. Integration & Verification

- [x] 4.1 Run targeted Clean 2D tests and relevant full test suite subset to verify no regressions after minimal adapters are connected.
- [x] 4.2 Add integration test running a complete clean-2d flow end-to-end, asserting final diagnostic state matches expected outcome for known molecules without changing expected geometry.
- [x] 4.3 Verify contract tests pass against actual module outputs (not just synthetic data).
- [x] 4.4 Verify normalized reporting score does not change existing ranking or candidate selection behavior.

## 5. Cleanup & Documentation

- [x] 5.1 Update README or docs to reference new diagnostic format and its usage.
- [x] 5.2 Remove only trivial local duplication created by this change; document non-trivial conceptual duplication as follow-up.
- [x] 5.3 Ensure internal metadata is preserved but ignored by contract tests where appropriate (per spec requirement).
