# Proposal: Repair local-template 3D publication and close viewport interaction gaps

## Problem

A user-local `Plantillas → Productos naturales → tetrandrina` template reaches RDKit with an aromatic bond to a terminal atom. RDKit rejects the graph during sanitization, before ETKDG, coordinate publication, scene construction, or OpenGL upload. The current build gesture previews a 2D screen-space line/circle rather than the 3D geometry that will be committed. Live MM frames can be delivered in a sub-frame burst and each callback requests a redraw, so Qt coalesces paints and the geometry appears to jump. The viewport also lacks a persistent camera-oriented axes indicator and explicit orbit guidance.

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

Validate and commit this OpenSpec first. Deliver the four implementation tracks as separate local commits, each with a failing regression test observed before production changes. Keep the OpenSpec active while real GUI acceptance remains pending.
