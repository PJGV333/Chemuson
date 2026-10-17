# Design

## Evidence and failure localization

The exact local template currently imports as a 43-atom, 48-bond MolGraph. It has six supported abbreviations (`CH3` ×2, `OCH3` ×4), and no explicit atom or bond stereochemistry. Expanding abbreviations on a copy yields 47 atoms and 52 bonds. The original graph has aromatic bridge bond 9 (atom 5–11) flagged aromatic while its stored order is single; atom 11 has degree one. `graph_conformer3d_stream` fails in RDKit sanitization with `non-ring atom 10 marked aromatic`; ETKDG is never called, so no coordinates reach scene construction. Clearing only this contradictory aromatic flag on a calculation copy, while retaining single-bond order and expanding the abbreviations, allows ETKDG to return coordinates. A prototype produced 93 scene atoms and a finite mesh, but MMFF did not converge; this is geometric availability, not a claim of minimized chemistry or visual acceptance.

## Architecture decisions

1. **Calculation-only molecular normalization.** Reuse `expand_abbreviations_for_calculation` on a defensive graph copy. Detect aromatic edges that are graph bridges. Only a bridge stored with an explicit single-bond order may be demoted from aromatic to single for RDKit calculation; preserve connectivity and all source IDs. Report every expansion/demotion in diagnostics. Leave unknown or ambiguous aromatic assignments as errors rather than guessing. Preserve all specified stereo fields and do not assign missing stereo.
2. **Derived scene atoms.** Keep the original canvas graph unchanged. Carry added abbreviation atoms/bonds as scene-only generated metadata, with deterministic IDs and coordinates, alongside RDKit-derived H. Coordinate sets keep calculation atom IDs so subsequent Live MM can optimize the same graph. Scene picking/editing must not route generated IDs into document mutation commands.
3. **3D DrawTool.** Track a world-space start/end point. Project the cursor onto the camera-facing plane through the source atom; when a real atom is hit, snap the endpoint to its exact 3D coordinates. Build a temporary render scene containing a tinted sphere and bond, separate from the pickable/document scene. On release, reuse existing atom/bond commands and valence checks; atom plus bond plus geometry belong to one undo macro. Escape removes preview without mutation.
4. **Double-click rule.** Double-clicking an atom changes that atom to the dock's selected build element, using the existing undoable element-change route and its valence/H safeguards. Double-clicking empty space extends the currently selected single-atom anchor by exactly one single bond; with no eligible anchor it does nothing and reports why. Delay only an empty-space single-click commit for Qt's double-click interval so the first click cannot create a duplicate before the second click is interpreted.
5. **Scientific frames versus visual playback.** Retain every RDKit `OptimizationFrame` as an immutable optimizer sample. Animate coordinates only between snapshots with identical atom IDs, force field, and generated-atom/bond identity. Interpolated coordinate sets carry no energy value and are tagged visual-only. `E` and `ΔE` remain measured snapshot values; `ΔE` is the difference between adjacent optimizer snapshots, never a gradient. The HUD shows force field, kcal/mol energy, sample ΔE, step/max steps, and the number of fixed atom-position constraints. Do not report an uncomputed gradient.
6. **Bounded playback and cancellation.** Use a bounded queue and a 30–60 Hz GUI timer; use a fixed visual transition interval when snapshots arrive in bursts. Editing, document switch, job replacement, or incompatible frame identity clears queued frames and requests worker cancellation before any later publication. Final coordinates remain the worker's scientific result; playback never writes interpolated coordinates to the document.
7. **Navigation overlays.** Draw the camera basis axes in a persistent lower-left viewport overlay. During empty-space orbit drags, draw transient directional/arc guidance and label camera motion; atom moves keep the existing geometry preview and receive a distinct status cue.

## Compatibility and license

No `.cmsn` serialization changes, new dependencies, or package ownership changes are planned. Avogadro source informs interaction and HUD behavior only; all implementation is written independently. No Avogadro source is copied.

## Commit boundaries

- Commit 1: this validated OpenSpec and baseline.
- Commit 2: exact local-template diagnostics/normalization and chemistry regression.
- Commit 3: world-space DrawTool preview, snapping, double-click, and undo tests.
- Commit 4: Live MM playback/HUD/cancellation and measured frame/animation tests.
- Commit 5: camera axes/orbit cues and navigation tests.

Each feature commit updates tasks and validation. Do not archive the change while GUI acceptance remains pending.
