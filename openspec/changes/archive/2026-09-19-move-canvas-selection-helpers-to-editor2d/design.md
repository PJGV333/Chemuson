# Design: Canonical editor2d selection ownership

## Package boundary

M20 owns the canonical implementations of `selection_geometry`,
`selection_bounds`, `selection_hit_testing`, `selection_overlay` and
`selection_clipboard` under `gui/editor2d/selection/`. The parent
`gui/editor2d/__init__.py` is a logic-free namespace file owned by M08. Future
siblings such as drawing, text, rendering and interactions are not M20-owned.

M09 remains the owner of canvas mixins and imports M20 for helper policy. The
five old `gui.canvas` module paths become thin star-import compatibility shims;
they contain no function or class definitions and preserve existing consumers.

## Catalog and documentation

`architecture/modules.yml` assigns only `gui/editor2d/selection/` and its tests
to M20, removes the five helpers from M09's internal API and test inventory,
and keeps M20 in M09's current and target dependencies. M20 has no temporary
exceptions or circular dependencies and never depends on M09.

`docs/modules/M09-canvas.md` documents the compatibility boundary and
`docs/modules/M20-editor2d-selection.md` documents canonical selection
ownership without claiming the whole editor2d namespace.

## Validation

Run the new architecture tests, all selection and clipboard regressions,
compileall, Ruff, OpenSpec strict validation, and the complete test suite.
Manual GUI smoke validation remains required before archiving.
