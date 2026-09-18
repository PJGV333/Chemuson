# Design: Canonical editor2d ownership

## Package boundary

M20 owns the canonical implementations of `selection_geometry`,
`selection_bounds`, `selection_hit_testing`, `selection_overlay` and
`selection_clipboard` under `gui/editor2d/`. The helpers retain their current
stdlib/PyQt6-only imports and therefore have no ChemUSON module dependencies.

M09 remains the owner of canvas mixins and imports M20 for helper policy. The
five old `gui.canvas` module paths become thin star-import compatibility shims;
they contain no function or class definitions and preserve existing consumers.

## Catalog and documentation

`architecture/modules.yml` assigns the editor2d package and its tests
exclusively to M20, removes the five helpers from M09's internal API and test
inventory, and adds M20 to M09's current and target dependencies. M20 has no
temporary exceptions or circular dependencies.

`docs/modules/M09-canvas.md` documents the compatibility boundary and
`docs/modules/M20-editor2d.md` documents canonical ownership and the stable
legacy imports.

## Validation

Run the new architecture tests, all selection and clipboard regressions,
compileall, Ruff, OpenSpec strict validation, and the complete test suite.
Manual GUI smoke validation remains required before archiving.
