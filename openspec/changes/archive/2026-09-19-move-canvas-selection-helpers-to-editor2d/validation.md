# Validation: Narrow M20 ownership to editor2d selection

## Automated validation

- M20 canonical paths and parent namespace ownership: 110 focused tests passed.
- `pytest -q tests/architecture/test_editor2d_selection_ownership.py
  tests/architecture/test_canvas_selection_geometry.py
  tests/architecture/test_canvas_selection_bounds.py
  tests/architecture/test_canvas_selection_hit_testing.py
  tests/architecture/test_canvas_selection_overlay.py
  tests/architecture/test_canvas_selection_clipboard.py`: 49 passed.
- `pytest -q tests/test_canvas_selection_geometry.py
  tests/test_canvas_selection_bounds.py tests/test_canvas_selection_hit_testing.py
  tests/test_canvas_selection_overlay.py
  tests/test_canvas_selection_clipboard_policy.py`: 65 passed.
- Architecture catalog/import/public API tests: 116 passed.
- Full architecture suite: 254 passed.
- Related clipboard, selection, diagram and canvas regressions: 112 passed.
- `python -m compileall -q src tests tools packaging`: passed.
- Ruff on migrated helpers, shims, consumers and architecture tests: passed.
- `openspec validate --all --strict`: 27 passed, 1 pre-existing unrelated
  failure in `spec/application-composition-root` because its first requirement
  lacks a SHALL or MUST keyword.
- `openspec validate move-canvas-selection-helpers-to-editor2d --strict`: valid.
- `git diff --check`: passed.

## Manual validation

Completed: Qt smoke testing exercised selection, overlay handles, copy/paste and
legacy import consumers. No behavioral regression was observed.
