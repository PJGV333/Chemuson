# Proposal: Repair local-template 3D publication and close viewport interaction gaps

## Problem

A user-local `Plantillas → Productos naturales → tetrandrina` template reaches RDKit with an aromatic bond to a terminal atom. RDKit rejects the graph during sanitization, before ETKDG, coordinate publication, scene construction, or OpenGL upload. A later identity audit traced this to the template itself: atom 11 duplicates atom 6's 2D position, is terminal, and is attached to atom 5 instead of the missing ring-closing edge 5–6. The import and `CH3`/`OCH3` expansion are faithful; the current aromatic-flag normalization permits a wrong `C39H46N2O6` graph (47 heavy atoms, 93 total) to reach a scene, rather than the reference `C38H42N2O6` graph (46 heavy, 88 total). Do not modify the user-local model without explicit authorization. The current build gesture previews a 2D screen-space line/circle rather than the 3D geometry that will be committed. Live MM frames can be delivered in a sub-frame burst and each callback requests a redraw, so Qt coalesces paints and the geometry appears to jump. The viewport also lacks a persistent camera-oriented axes indicator and explicit orbit guidance.

## Change

Preserve the authored 2D MolGraph while using a diagnosed calculation copy for 3D conformers; render any expanded abbreviation atoms as derived scene geometry. Implement a true world-space draw preview with target snapping and atomic undo, define contextual double-click behavior, animate only the visual path between compatible optimizer snapshots, expose scientifically labeled telemetry, and add camera axes/orbit cues.

## Scope and safeguards

- Keep all work in existing ChemUSON packages and `feature/chemuson-3d`.
- Do not infer or add stereochemistry absent from the local template.
- Do not mutate the saved template or source MolGraph during calculation normalization.
- Preserve valence checks, explicit/derived H behavior, undo/redo, and current dock styling.
- Use Avogadro 1.2 source only as behavioral evidence; do not copy GPLv2-only code.
- Do not claim a visible render from `render_requested`, offscreen widgets, or mocks. If a real GUI session is not authorized, leave the manual guide and report the visual limitation.
- No push, merge, or PR.

## Delivery

Validate and commit the OpenSpec before production changes. Keep the prior specification, chemistry-normalization, and direct-drawing work in their separate local commits. Group the remaining Live MM and navigation UI integration in one descriptive viewport commit because both use the same widget, then commit the identity audit and manual acceptance record separately. Keep this OpenSpec active while the source-template identity and real GUI/OpenGL acceptance remain pending; do not push, merge, or open a PR.
