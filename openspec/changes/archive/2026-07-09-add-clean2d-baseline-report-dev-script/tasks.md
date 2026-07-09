# Tasks: Add Clean 2D Baseline Report Developer Script

## 1. OpenSpec Artifacts
- [x] Create proposal, design, tasks, and spec delta.
- [x] Validate `add-clean2d-baseline-report-dev-script` before implementation.

## 2. Script Implementation
- [x] Add `scripts/clean2d_baseline_report.py`.
- [x] Implement `write` command.
- [x] Implement `compare` command.
- [x] Implement `review` command.
- [x] Return explicit exit codes for equivalent, changed, and error outcomes.
- [x] Keep outputs JSON-stable and observational-only.

## 3. Tests
- [x] Verify `write` creates a valid JSON report.
- [x] Verify compare equivalent reports exits `0`.
- [x] Verify compare changed reports exits `1`.
- [x] Verify review equivalent reports exits `0`.
- [x] Verify review changed reports exits `1` and prints review JSON.
- [x] Verify missing or invalid JSON inputs exit `2`.
- [x] Verify the script remains developer/test-only and does not activate gates.

## 4. Validation And Commit
- [x] Run requested OpenSpec validations and pytest suites.
- [x] Commit as `test(clean2d): add baseline report developer script`.
