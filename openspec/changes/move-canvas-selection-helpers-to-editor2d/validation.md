# Validation: Move canvas selection helpers to editor2d

## Automated validation

- Red architecture test before migration: 5 failures, caused by absent M20
  canonical paths, non-shim legacy modules and missing M20 catalog entry.
- `pytest -q tests/architecture/test_editor2d_selection_ownership.py
  tests/architecture/test_canvas_selection_geometry.py
  tests/architecture/test_canvas_selection_bounds.py
  tests/architecture/test_canvas_selection_hit_testing.py
  tests/architecture/test_canvas_selection_overlay.py
  tests/architecture/test_canvas_selection_clipboard.py`: 42 passed.
- `pytest -q tests/test_canvas_selection_geometry.py
  tests/test_canvas_selection_bounds.py tests/test_canvas_selection_hit_testing.py
  tests/test_canvas_selection_overlay.py
  tests/test_canvas_selection_clipboard_policy.py`: 65 passed.
- Architecture catalog/import/public API tests: 116 passed.
- Related clipboard, selection, diagram and canvas regressions: 112 passed.
- `python -m compileall -q src tests tools packaging`: passed.
- Ruff on migrated helpers, shims, consumers and architecture tests: passed.
- `openspec validate move-canvas-selection-helpers-to-editor2d --strict`: valid.
- `git diff --check`: passed.

## Manual validation

Pending: launch the Qt application and exercise selection, overlay handles,
copy/paste and legacy import consumers before archiving this OpenSpec.
