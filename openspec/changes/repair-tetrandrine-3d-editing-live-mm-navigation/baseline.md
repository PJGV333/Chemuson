# Baseline

## Repository state
- Branch: `feature/chemuson-3d`
- HEAD: `805a09a543eba801a662e79f748267fe4f58ed19` (`feat(geometry3d): harden complex conformer generation`)
- Worktree: clean before this OpenSpec was created.
- Existing floating-dock fix is preserved at `042eed4`; no remote operations were performed.

## Required baseline checks
- `./.venv/bin/python -m compileall -q src tests tools packaging` — exit 0, no output.
- `./.venv/bin/pytest --collect-only -q` — exit 0, 2,243 tests collected.
- Focused baseline:
  `./.venv/bin/pytest -q tests/test_viewport3d_widget.py tests/test_viewport3d_camera.py tests/test_viewport3d_geometry.py tests/test_geometry3d_service.py tests/test_rdkit_live_stream.py tests/test_compchem3d_dock.py tests/test_template_library.py`
  — exit 0, 164 passed, 2 skipped in 48.83 s.
- `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741` — exit 1 for the existing `F401` unused `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Full `pytest -q` was not rerun: the repository's known Qt/SIP SIGSEGV remains open. Collection and focused baseline were run instead.

## Exact local tetrandrina template reproduction
The app's user-local library was read without instantiating the write-on-load `TemplateLibrary` and without changing its file. The selected entry is `Productos naturales / tetrandrina`.

- Imported exact MolGraph: 43 atoms, 48 bonds; element/alias counts: C 33, CH3 2, N 2, O 2, OCH3 4.
- Explicit atom stereochemistry: none. Explicit bond stereochemistry: none. No stereochemistry was inferred.
- `expand_abbreviations_for_calculation` on a copy: 47 atoms, 52 bonds; C 39, N 2, O 6.
- Original 3D worker failure: `invalid_graph`, stage `sanitization`, detail `sanitize_failed:non-ring atom 10 marked aromatic`. ETKDG and coordinate publication were not reached.
- Graph diagnosis: aromatic bridge bond ID 9 joins atom IDs 5 and 11; atom 11 has degree 1; stored bond order is single. The flag is internally contradictory because an acyclic bridge cannot be aromatic.
- Read-only prototype on an expanded copy, demoting only that bridge's aromatic flag while retaining its single order: ETKDG return code 0; embedding 0.218 s; minimization 0.097 s but did not converge; result still carried coordinates. Scene construction yielded 93 atoms and a 120,312-vertex finite mesh. This verifies the chemistry-to-mesh path, not OpenGL visibility.

## Follow-up baseline before the explicitly authorized local-template repair (2026-10-17)

- Branch/HEAD: `feature/chemuson-3d`, `051936c757c311ab369e0c9796d837ea4ec64702`; worktree was clean before this OpenSpec scope update.
- User-library file: `~/.config/Chemuson/template_library.json`, 10,844 bytes; SHA-256 before repair: `75051c9bc73e92d34bc53f109a123aa909cbf1e47186f8e9f6dd7fafdab5c68d`.
- Read-only selection found exactly one `Productos naturales / tetrandrina` entry, ID `tpl_f79a6425241340fab69f51e0b9d65870`. No ChemUSON process was running.
- `./.venv/bin/python -m compileall -q src tests tools packaging` — exit 0, no output.
- `./.venv/bin/python -m pytest --collect-only -q` — exit 0, 2,265 tests collected in 0.47 s.
- Focused existing chemistry baseline `./.venv/bin/python -m pytest -q tests/test_geometry3d_aromaticity.py tests/test_geometry3d_service.py tests/test_clean2d_tetrandrina_like_local.py` — 44 passed, 1 skipped in 6.47 s.
- `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741` — exit 1 for the pre-existing unused `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Full pytest was not run because the repository's documented Qt/SIP SIGSEGV remains open.
- In-memory chemistry prototype only: the proposed 5–11 → 5–6 correction sanitized to `C38H42N2O6`, 46 heavy atoms, 52 heavy bonds, 1 component, cycle rank 7, 24 aromatic atoms/bonds; its InChI matched the PubChem CID 73078 connectivity with stereochemistry unspecified. No config or tracked production file was changed during this prototype.

## Live MM worker-frame baseline
Measured with real RDKit MMFF94 streamed frames from collinear starting coordinates; timestamps are Python consumer receipt times, not display or paint timestamps.

- Ethane, 250 max steps / 25 steps per update: 2 samples; total 0.124690 s; first sample 0.110921 s; the two energies were both -4.734365 kcal/mol.
- Butane, same settings: 8 samples over 0.139618 s; first sample 0.122623 s; seven inter-sample intervals ranged 0.327–0.456 ms (median 0.379 ms), all below a 16 ms frame.
- Butane's measured E sequence (kcal/mol): step 25 `89.368972`, 50 `-2.841572`, 75 `-4.239516`, 100 `-4.292270`, 125 `-4.293689`, 150 `-4.293761`, 175 `-4.293764`, 200 `-4.293764`.
- Current `_on_compchem_frame_ready` refreshes the viewport immediately for each callback. Eight samples arriving within a few milliseconds are faster than a display frame, and queued `update()` paints can coalesce; this explains the apparent jump. It is not evidence that RDKit emitted no intermediate samples.

## Avogadro 1.2 source inspection
Read-only inspection of `DrawTool` and `AutoOptTool` in the local Avogadro 1.2 source checkout. `DrawTool` tracks a temporary end atom/bond, snaps to existing hit atoms, positions endpoints via camera unprojection, and groups related edits on release. `AutoOptTool` paints E/ΔE and constraint count; its 50 ms timer calls `AutoOptThread::update()` directly, whose implementation performs force-field setup/steps in that call path. The source is GPLv2; only behavior was studied and no code copied. The local source is sufficient evidence; REA/Ghidra was not needed for these source-level questions.

## Baseline conclusion
The local template's first blocking failure is RDKit sanitization, before ETKDG/scene/OpenGL. The Live MM discontinuity is explained by measured sub-frame frame arrivals plus immediate scene refresh. Visual OpenGL acceptance remains untested; `render_requested` is not considered proof of visible output.

## Post-repair verification (2026-10-17)

- Local library entry ID `tpl_f79a6425241340fab69f51e0b9d65870` was changed only in its `molblock` and `smiles` values. Full backup `~/.config/Chemuson/template_library.json.bak-20261017T193637027947-0700` was verified byte-for-byte (10,844 bytes; SHA-256 `75051c9bc73e92d34bc53f109a123aa909cbf1e47186f8e9f6dd7fafdab5c68d`). Updated JSON is 10,778 bytes; SHA-256 `0f34a831b22577c1d954788638cf6950dd53cb2960cbd2c533b31446c474eb14`. All other serialized bytes and parsed library fields/entries were verified unchanged; neither config nor backup is in Git.
- Target MolGraph: 42 atoms / 48 bonds; expanded copy: 46 heavy atoms / 52 bonds, `C38H42N2O6`, one component, cycle rank 7, 24 aromatic atoms/bonds. InChI matched the CID 73078 connectivity; regenerated stored SMILES round-tripped to the same InChI. No stereo was added.
- Corrected fixture was derived by a targeted MolBlock record edit: delete atom 11, replace 5–11 with aromatic 5–6, adjust remaining atom-index references, preserve other coordinates/attributes. Historical fixture remains unchanged.
- Worker probe (seed 173, MMFF94, max_iters 80): stage `done`, 46 finite coordinates; minimization did not converge (`method=unoptimized`, `energy=None`), total 0.108235 s. CPU scene: 88 atoms / 94 bonds and finite mesh; this is not OpenGL acceptance.
- Post-repair focused suite: 48 passed, 1 skipped in 7.55 s. Compileall and scoped Ruff on changed Python files passed. Broad selected Ruff still reports the pre-existing unused `math` import at `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. Full pytest was not run because the known Qt/SIP SIGSEGV remains open.
- Local code commit: `be7eca1` (`fix(geometry3d): reject invalid aromatic bridges`). No push, merge, or PR. Actual GUI/OpenGL acceptance remains pending.
