# Manual GUI/OpenGL acceptance guide

Use this only in an authorized desktop session with the project environment. Do not count an offscreen/mock test as visual acceptance. Capture the application version, OpenGL vendor/renderer, viewport status and screenshots for each result.

## 1. Local tetrandrina template
1. Open the actual local ChemUSON template library and insert `Plantillas → Productos naturales → tetrandrina`.
2. Open/show the 3D dock and wait for conformer generation.
3. Confirm atoms and bonds are visibly rendered (not a black/empty viewport); rotate and fit the molecule.
4. Confirm the status/diagnostics report abbreviation expansion, the aromatic bridge normalization, and any MMFF nonconvergence without inventing stereo.
5. Confirm the authored 2D template remains unchanged and generated scene atoms are not treated as document atoms.

## 2. DrawTool and edit history
1. Enable 3D editing and choose `Construir molécula`.
2. Drag from an atom into empty space; verify a world-space sphere and bond follow the cursor. Release and verify one new atom and one bond.
3. Undo once: both atom and bond disappear. Redo once: both return.
4. Drag from one atom to another; verify visual snap and exactly one bond, with no duplicate endpoint atom.
5. Double-click an atom: selected element replaces it; verify valence/H behavior and one undoable action.
6. Select one atom and double-click empty space: exactly one single-bonded chain atom appears. Verify the first click does not create an extra atom.
7. Start a drag and press Escape; verify no graph, coordinate, or undo-stack mutation.
8. Repeat representative edit, H, valence, undo/redo cases with ethane and butane.

## 3. Live MM
1. Start RDKit Live MM for ethane, then butane, from a non-minimized 3D arrangement.
2. Confirm multiple geometric transitions are perceptible even when worker samples arrive together.
3. Confirm HUD fields: force field, measured `E` and sample-to-sample `ΔE` in kcal/mol, current/max step, and fixed atom-position count. Do not expect a gradient field unless a separately computed gradient is explicitly labeled.
4. Edit the molecule during optimization; verify the displayed queue clears immediately, the stale worker is cancelled, and no old frame repopulates the scene.
5. Confirm interpolation is visual-only and the saved/final coordinate set equals an optimizer snapshot, never an interpolated point.

## 4. Navigation
1. Rotate the camera around ethane and butane; verify persistent lower-left X/Y/Z axes rotate with the view.
2. Drag in empty space; verify a transient orbit cue and camera-motion status.
3. Move an atom in edit mode; verify a distinct atom-motion cue and no orbit indicator.
4. Confirm overlays do not affect picking, selection, chemistry, or geometry.
