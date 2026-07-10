## 1. Analysis and Preparation

- [ ] 1.1 Analyze `src/chemuson/clean2d/engine.py` to identify the `smart propose` flow and where the repair check should be integrated.
- [ ] 1.2 Identify the existing repair engine entry points used by the system.

## 2. Core Implementation

- [ ] 2.1 Modify `smart propose` flow to check for distortion and apply a safe repair before proposing alternatives.
- [ ] 2.2 Ensure the repaired geometry is passed correctly to the alternative generation logic.

## 3. Testing and Verification

- [ ] 3.1 Update `tests/test_clean2d_smart_propose.py` to verify that the repair is applied correctly.
- [ ] 3.2 Verify that the minimum bond length requirements are met in the proposed alternatives.
- [ ] 3.3 Run existing tests to ensure no regressions in the `clean2d` engine.
