# Proposal: Narrow M20 ownership to editor2d selection

## Why

The extracted selection helpers are deterministic editor-2D policy, but M20
must own only the selection namespace rather than the entire future `gui.editor2d`
package. The parent namespace must remain available for future sibling modules.

## Scope

- Add `src/chemuson/gui/editor2d/selection/` as the canonical M20 package.
- Move the five selection helper implementations below that namespace.
- Keep `src/chemuson/gui/editor2d/__init__.py` as a logic-free parent namespace
  owned by M08.
- Leave compatibility shims at the old `gui.canvas` paths.
- Update imports, tests, module catalog and documentation.

## Non-goals

- Do not change helper behavior or payload formats.
- Do not move canvas mixins, Qt event dispatch, scene mutation, commands or undo.
