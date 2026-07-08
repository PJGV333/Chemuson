# Design: Clean 2D Baseline Report Developer Script

## Overview
The script is a developer/test-only wrapper around existing test-owned baseline report helpers. It lives under `scripts/` and imports only the test helper modules that already generate, compare, and review Clean 2D baseline reports.

## Commands
### write
`python scripts/clean2d_baseline_report.py write --output report.json`

Generates a full Clean 2D baseline report using `build_baseline_report`, serializes it with `baseline_report_to_json`, creates parent directories for the requested path, and writes only to that explicit path.

### compare
`python scripts/clean2d_baseline_report.py compare --left old.json --right new.json`

Loads two report JSON files, runs `compare_baseline_reports`, prints stable JSON diff output, and exits `0` when equivalent or `1` when different.

### review
`python scripts/clean2d_baseline_report.py review --left old.json --right new.json`

Loads two report JSON files, compares them, classifies the diff with `classify_baseline_diff`, summarizes it with `summarize_baseline_diff_review`, prints stable JSON summary output, and exits `0` when no review items exist or `1` when items exist.

## Exit Codes
- `0`: successful write, equivalent compare, or no review items.
- `1`: observable differences or review items exist.
- `2`: argument, read, write, or JSON error.

Exit code `1` is observational in this change; it is not a geometry gate.

## Test Boundary
The script does not modify production code or public document formats. It is intentionally a local developer utility for test-owned baseline artifacts.

## Error Handling
The script reports errors to stderr and returns `2` for invalid arguments, missing files, invalid JSON, or write failures.
