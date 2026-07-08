## 1. Contract Surface

- [x] 1.1 Review current Clean 2D result objects, controller attempts, candidate metadata, and status messages against the spec vocabulary.
- [x] 1.2 Identify the smallest reporting surface that can expose `applied`, `rejected`, `preserve-only`, `no-op`, and `failed-controlled` without changing layout algorithms.
- [x] 1.3 Map existing internal rejection messages and safety failures to the stable reasons `invalid-coordinates`, `invariant-violation`, `stereo-risk`, `boundary-bond-risk`, `new-crossing-risk`, `collision-risk`, `collapsed-ring-risk`, `excessive-displacement`, `worse-quality`, and `backend-failure`.
- [x] 1.4 Document any current behavior that cannot yet be mapped cleanly and keep it as diagnostic detail rather than changing algorithms.

## 2. Contract Tests

- [x] 2.1 Add or update tests proving Clean 2D preserves atom IDs, bond IDs, atom count, bond count, connectivity, bond orders, charges, labels, stereo metadata, and selection.
- [x] 2.2 Add or update tests proving rejected or unavailable candidates leave coordinates and non-coordinate data unchanged.
- [x] 2.3 Add or update tests proving selection-targeted Clean 2D validates boundary bonds against the full canvas model.
- [x] 2.4 Add or update tests proving applied coordinate changes are a single undoable operation and rejected candidates do not push coordinate-changing undo commands.

## 3. Safety Reason Tests

- [x] 3.1 Add tests for stable `invalid-coordinates` reporting when candidate coordinates are missing or non-finite.
- [x] 3.2 Add tests for stable `invariant-violation` reporting when a candidate changes non-coordinate chemical data.
- [x] 3.3 Add tests for stable `stereo-risk` and `boundary-bond-risk` reporting where current validation support exists.
- [x] 3.4 Add tests for stable `new-crossing-risk`, `collision-risk`, `collapsed-ring-risk`, `excessive-displacement`, and `worse-quality` reporting where current safety support exists.
- [x] 3.5 Add tests or controlled mocks for stable `backend-failure` reporting when an optional backend fails without mutating the structure.

## 4. Mode and Candidate Diagnostics

- [x] 4.1 Add tests that `quick`, `publication`, and `propose` retain documented intent at the observable level without asserting exact backend order.
- [x] 4.2 Add tests that evaluated candidates expose reportable source labels for diagnostics.
- [x] 4.3 Ensure backend-specific detail can remain diagnostic metadata while stable result state and reason values remain testable.

## 5. Validation

- [x] 5.1 Run the focused Clean 2D, selection, and geometry test suites affected by the contract work.
- [x] 5.2 Run OpenSpec validation for `document-clean2d-decision-contract`.
- [x] 5.3 Record any implementation gaps found during tests as follow-up OpenSpec changes rather than expanding this change into algorithm work.

## Implementation Notes

- Stable reporting was added through existing result/attempt/candidate objects and reason mapping helpers; no layout algorithm, MolGraph, canvas, or document format changes were introduced.
- Focused contract/reporting tests pass. Running the broader `tests/test_clean2d_safety.py` suite revealed `test_run_clean_2d_cyclic_geometry_polish_repairs_double_bond_linker` failing in isolation with a double-bond linker length mismatch. This is a geometry behavior gap and remains out of scope for this contract-only change.
