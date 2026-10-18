# Manual GUI/OpenGL acceptance guide

## CachyOS launch

Run this only in an authorized, logged-in desktop session; no real GUI/OpenGL acceptance was performed for this campaign:

```sh
cd /home/ccachyavgp/Documentos/Chemuson-3D
./.venv/bin/python -m chemuson
```

Use the actual desktop viewport and hardware OpenGL path. Record the app version, OpenGL vendor/renderer, viewport status and screenshots. Do not count mocks, offscreen tests, scene meshes, `render_requested`, or upload success as visible acceptance.

## 1. Local tetrandrina identity gate

The local entry has been repaired and read back from `~/.config/Chemuson/template_library.json`: its expanded identity is `C38H42N2O6` (46 heavy / 88 total), with one connected component and connectivity matching PubChem CID 73078 without stereochemistry. The canonical fallback SMILES was regenerated and round-tripped to the same connectivity. A full backup was verified before the atomic update; only this entry's MolBlock and SMILES fields changed.

The worker produced 46 finite heavy-atom coordinates, but MMFF94 did not converge (`method=unoptimized`, `energy=None`). Treat the geometry as an unoptimized conformer, not a minimized structure. CPU scene/mesh tests are not visual acceptance.

In a fresh, authorized CachyOS desktop session, reopen `Plantillas → Productos naturales → tetrandrina` and confirm the loaded structure, formula, components, connectivity, valences, and absence of invented stereo. Inspect the actual hardware-OpenGL viewport for visible atoms/bonds, fit and rotate it, and verify generated scene atoms do not mutate the authored 2D graph. Record application version, OpenGL vendor/renderer, viewport status, and screenshots. This campaign has **not** verified OpenGL pixels; acceptance remains pending.

## 2. DrawTool and edit history
1. Open ethane or butane, show the 3D dock, enable 3D editing, and choose `Construir molécula`.
2. Drag from an atom into empty space; verify a world-space sphere and bond follow the cursor. Release and verify one new atom and one bond.
3. Undo once: both atom and bond disappear. Redo once: both return.
4. Drag from one atom to another; verify visual snap and exactly one bond, with no duplicate endpoint atom.
5. Double-click an atom: selected element replaces it; verify valence/H behavior and one undoable action.
6. Select one atom and double-click empty space: exactly one single-bonded chain atom appears. Verify the first click does not create an extra atom.
7. Start a drag and press Escape; verify no graph, coordinate, or undo-stack mutation.

## 3. Live MM
1. Start RDKit Live MM for ethane, then butane, from a non-minimized 3D arrangement.
2. Confirm multiple geometric transitions are perceptible even when worker samples arrive together and are presented in step order.
3. Confirm HUD fields: force field, measured `E` and sample-to-sample `ΔE` in kcal/mol, current/max step, and fixed atom-position count. Energy values must change only with real samples; interpolated coordinates are labeled visual-only and carry no scientific energy/gradient.
4. Edit the molecule during optimization; verify the displayed queue clears immediately, the stale worker is cancelled, and no old frame repopulates the scene.
5. After playback drains, verify the final document coordinates equal the worker's final optimizer snapshot, never an interpolated point.

## 4. Navigation
1. Rotate the camera around ethane and butane; verify persistent lower-left X/Y/Z axes rotate with the view.
2. Drag in empty space; verify a transient orbit cue and camera-motion status.
3. Move an atom in edit mode; verify a distinct atom-motion cue and no orbit indicator.
4. Confirm overlays do not affect picking, selection, chemistry, or geometry.
