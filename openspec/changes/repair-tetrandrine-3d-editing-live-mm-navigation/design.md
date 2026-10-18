# Design

## Evidence and failure localization

Before the authorized repair, the local template imported as a 43-atom, 48-bond MolGraph with six supported abbreviations (`CH3` ×2, `OCH3` ×4) and no explicit atom or bond stereochemistry. Expanding abbreviations on a copy yielded 47 atoms and 52 bonds. The original graph had aromatic bridge bond 9 (atom 5–11) flagged aromatic while its stored order was single; atom 11 had degree one. `graph_conformer3d_stream` failed during RDKit sanitization with `non-ring atom 10 marked aromatic`; ETKDG was never called. An earlier implementation cleared that flag on a calculation copy and produced coordinates, but identity auditing proved this was an invalid workaround: it represented the wrong molecule. The normalizer has been replaced with an explicit pre-worker rejection for acyclic aromatic bridges. The historical prototype's finite mesh and nonconverged MMFF result did not establish molecular identity or visual acceptance.

## Post-implementation identity audit

A deeper audit shows that the sanitization failure is a symptom of a malformed user-local graph, not an isolated RDKit aromaticity issue. The MolBlock and its fixture are identical, and import/abbreviation expansion preserve the source connectivity. Source atom 11 duplicates atom 6's 2D coordinates, is terminal, and is connected to atom 5 by the invalid aromatic bridge; the expected ring-closing edge 5–6 is absent. Under the earlier flag normalizer, RDKit treated atom 11 as a methyl group, yielding `C39H46N2O6` (47 heavy atoms, 93 total), not the reference `C38H42N2O6` (46 heavy atoms, 88 total). That normalized graph passed ordinary valence checks, proving sanitization alone is not a molecular-identity check. The corrected identity and preserved pre-repair evidence are recorded in `identity-audit.md`.

The owner authorized repairing only the saved local `Productos naturales / tetrandrina` entry. The change is complete: atom 11 and its 5–11 bond were removed, aromatic closure 5–6 restored, and the stored SMILES regenerated from the expanded corrected graph. The original malformed fixture remains as historical evidence. A verified full backup was made and the JSON was atomically replaced by changing only the target MolBlock and SMILES string values; all other bytes and parsed library data were preserved. The corrected entry matches PubChem CID 73078 connectivity without assigning absent stereochemistry.

## Architecture decisions

1. **Strict aromatic bridge handling.** Expand abbreviations on a defensive calculation copy, but never demote an acyclic aromatic flag to make an invalid graph sanitize. Return a stage-specific diagnostic naming aromatic bridge bond IDs and preserve the source graph. Valid aromatic rings remain unchanged. Preserve all specified stereo fields and do not assign missing stereo.
2. **Derived scene atoms.** Keep the original canvas graph unchanged. Carry added abbreviation atoms/bonds as scene-only generated metadata, with deterministic IDs and coordinates, alongside RDKit-derived H. Coordinate sets keep calculation atom IDs so subsequent Live MM can optimize the same graph. Scene picking/editing must not route generated IDs into document mutation commands.
3. **3D DrawTool.** Track a world-space start/end point. Project the cursor onto the camera-facing plane through the source atom; when a real atom is hit, snap the endpoint to its exact 3D coordinates. Build a temporary render scene containing a tinted sphere and bond, separate from the pickable/document scene. On release, reuse existing atom/bond commands and valence checks; atom plus bond plus geometry belong to one undo macro. Escape removes preview without mutation.
4. **Double-click rule.** Double-clicking an atom changes that atom to the dock's selected build element, using the existing undoable element-change route and its valence/H safeguards. Double-clicking empty space extends the currently selected single-atom anchor by exactly one single bond; with no eligible anchor it does nothing and reports why. Delay only an empty-space single-click commit for Qt's double-click interval so the first click cannot create a duplicate before the second click is interpreted.
5. **Scientific frames versus visual playback.** Retain every RDKit `OptimizationFrame` as an immutable optimizer sample. Animate coordinates only between snapshots with identical atom IDs, force field, and generated-atom/bond identity. Interpolated coordinate sets carry no energy value and are tagged visual-only. `E` and `ΔE` remain measured snapshot values; `ΔE` is the difference between adjacent optimizer snapshots, never a gradient. The HUD shows force field, kcal/mol energy, sample ΔE, step/max steps, and the number of fixed atom-position constraints. Do not report an uncomputed gradient.
6. **Bounded playback and cancellation.** Use a bounded queue and a 30–60 Hz GUI timer; use a fixed visual transition interval when snapshots arrive in bursts. Editing, document switch, job replacement, or incompatible frame identity clears queued frames and requests worker cancellation before any later publication. Final coordinates remain the worker's scientific result; playback never writes interpolated coordinates to the document.
7. **Navigation overlays.** Draw the camera basis axes in a persistent lower-left viewport overlay. During empty-space orbit drags, draw transient directional/arc guidance and label camera motion; atom moves keep the existing geometry preview and receive a distinct status cue.

## Compatibility and license

No `.cmsn` serialization changes, new dependencies, or package ownership changes are planned. Avogadro source informs interaction and HUD behavior only; all implementation is written independently. No Avogadro source is copied.

## Commit boundaries

- Commit 1: validated OpenSpec and baseline (already committed).
- Commit 2: calculation-copy aromaticity normalization and original-graph regression (already committed); this fixes RDKit sanitization only, not the local template's chemical identity.
- Commit 3: world-space DrawTool preview, snapping, double-click, and undo tests (already committed).
- Commit 4: remaining viewport integration—timed Live MM playback/HUD/cancellation plus camera axes/orbit cues and interaction tests, grouped because these changes share the viewport widget.
- Commit 5: prior identity audit, CachyOS manual acceptance guide, and validation record.
- Commit 6: completed as `be7eca1` — authorized local-repair fixture, invalid-aromatic-bridge refusal, and focused regression tests. The user configuration and backup remain outside Git.
- Commit 7: final identity audit and validation evidence for the repaired entry, recorded separately from the code commit while preserving original defect history.

Do not archive the change while real GUI/OpenGL acceptance remains pending.
