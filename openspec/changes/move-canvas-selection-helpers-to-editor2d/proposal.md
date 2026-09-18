# Proposal: Move canvas selection helpers to editor2d

## Why

The extracted selection helpers are deterministic editor-2D policy, but they remain
owned by the M09 canvas package. A dedicated M20 package makes the ownership
boundary explicit while preserving the historical canvas import paths.

## Scope

- Add `src/chemuson/gui/editor2d/` as module M20.
- Move the five selection helper implementations there.
- Leave compatibility shims at the old `gui.canvas` paths.
- Update imports, tests, module catalog and documentation.

## Non-goals

- Do not change helper behavior or payload formats.
- Do not move canvas mixins, Qt event dispatch, scene mutation, commands or undo.
