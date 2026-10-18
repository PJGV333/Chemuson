# Proposal: Repair local-template 3D publication and close viewport interaction gaps

## Problem

At campaign start, the user-local `Plantillas → Productos naturales → tetrandrina` template reached RDKit with an aromatic bond to a terminal atom. RDKit rejected the graph during sanitization, before ETKDG, coordinate publication, scene construction, or OpenGL upload. The identity audit traced this to the template itself: atom 11 duplicated atom 6's 2D position, was terminal, and was attached to atom 5 instead of the missing ring-closing edge 5–6. Import and `CH3`/`OCH3` expansion were faithful; the earlier aromatic-flag normalizer let the wrong `C39H46N2O6` graph (47 heavy atoms, 93 total) reach a scene instead of the reference `C38H42N2O6` graph (46 heavy, 88 total). After explicit owner authorization, the local entry has now been narrowly repaired by removing atom 11 and replacing the invalid 5–11 bridge with aromatic closure 5–6; all other library data was preserved. The campaign also addressed the build gesture's 2D screen-space preview, sub-frame Live MM bursts, and missing camera-oriented axes/orbit cues. Real GUI/OpenGL acceptance remains pending.

## Change

Repair only the authorized local tetrandrina entry, preserving every other atom, substituent, bond, coordinate, category, ID, preference, and template. Keep the authored graph unchanged during conformer calculation; reject invalid aromatic bridge flags rather than silently de-aromatizing them. Render expanded abbreviation atoms as derived scene geometry. Implement a true world-space draw preview with target snapping and atomic undo, define contextual double-click behavior, animate only the visual path between compatible optimizer snapshots, expose scientifically labeled telemetry, and add camera axes/orbit cues.

## Scope and safeguards

- Keep all work in existing ChemUSON packages and `feature/chemuson-3d`.
- Do not infer or add stereochemistry absent from the local template.
- Do not mutate the saved template or source MolGraph during conformer calculation.
- Modify only the authorized tetrandrina MolBlock and its derived SMILES in the user JSON; preserve all other bytes where practical, use a verified full backup and atomic replacement, and never instantiate the write-on-load library for editing.
- Preserve valence checks, explicit/derived H behavior, undo/redo, and current dock styling.
- Use Avogadro 1.2 source only as behavioral evidence; do not copy GPLv2-only code.
- Do not claim a visible render from `render_requested`, offscreen widgets, or mocks. If a real GUI session is not authorized, leave the manual guide and report the visual limitation.
- No push, merge, or PR.

## Hito A — OpenGL fault isolation and repair

The owner reports that the corrected tetrandrina appears in the 2D canvas while the 3D viewport is black except for the XYZ overlay and `render_requested`. Do not change the repaired graph or template library. Trace the complete path from `main_window.py` publication through `SceneMolecule3D`, CPU mesh, VBO/VAO, shaders, camera matrices/clipping, draw submission, Qt framebuffer, and dock context lifecycle. Add stage-specific diagnostics and regression tests before changing the rendering path. Verify a synthetic molecule and the exact corrected tetrandrina through a real KDE Qt/OpenGL framebuffer readback; a mesh, upload flag, or draw submission alone is not visual acceptance. Preserve a truthful distinction between draw submission and non-background molecular pixels.

## Delivery

Validate and commit the OpenSpec before the authorized chemistry repair. Keep prior specification, conformer, direct-drawing, and Live MM/navigation work in separate local commits. Add the local template correction with focused regression tests and a separate descriptive commit; commit the identity audit/manual acceptance/validation record separately. Complete Hito A before writing Hito B or Hito C architecture changes. Keep this OpenSpec active until real GUI/OpenGL acceptance is evidenced; do not push, merge, or open a PR.
