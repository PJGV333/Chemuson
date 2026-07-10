# Tasks: Fix smart propose repair

- [ ] **Investigation**: Verify the exact conditions in `rank_clean2d_candidates` that trigger the rejection for distorted geometries.
- [ ] **Investigation**: Identify the best "safe repair" candidate to use when the geometry is bad.
- [ ] **Implementation**: Modify `rank_clean2d_candidates` to apply the repair before the rejection logic.
- [ ] **Implementation**: Ensure the repair doesn't violate any chemical invariants.
- [ ] **Testing**: Update `tests/test_clean2d_smart_propose.py` to verify the repair flow.
- [ ] **Validation**: Run the full suite and OpenSpec validation.
