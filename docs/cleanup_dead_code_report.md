# Cleanup dead code report

Branch: `cleanup/remove-dead-code-and-tests`

## Summary

Poda conservadora centrada en basura histórica y bajo riesgo: scripts manuales fuera del producto, una fachada legacy muerta, imports muertos, hacks `sys.path` en tests y algunos tests PR renombrados a suites semánticas. No se tocó lógica Clean2D/RDKit/ChemName crítica salvo imports muertos y nombres de tests.

## Baseline before cleanup

| Command | Result | Notes |
| --- | --- | --- |
| `python -m compileall src` | Passed | Source tree compiled. |
| `pytest -q` | Failed | 842 passed, 55 skipped, 6 failed in ~53s. Failures were already present on `main`. |
| `pytest --collect-only -q` | Passed | 903 tests collected in ~2.4s. |

Baseline failing tests:

| Test | Area | Observed failure |
| --- | --- | --- |
| `tests/test_clean2d_aromatic_mixed.py::test_naphthalene_fused_system_does_not_collapse` | Clean2D | Engine returns `complex_preserve_unsafe` instead of `ok`. |
| `tests/test_clean2d_multilayer_constraints.py::test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph` | Clean2D | No selected candidate; engine preserves complex projection. |
| `tests/test_clean2d_safety.py::test_run_clean_2d_cyclic_geometry_polish_repairs_double_bond_linker` | Clean2D | Double bond linker length remains too long. |
| `tests/test_clean2d_smart_propose.py::test_smart_clean_repairs_distorted_structure_before_proposing` | Clean2D | Minimum post-clean bond length is below expected threshold. |
| `tests/test_clean2d_terminal_rings.py::test_single_anchor_logic_does_not_run_for_fused_rings` | Clean2D | Engine returns `complex_preserve_unsafe` instead of `ok`. |
| `tests/test_clean2d_tetrandrine_like_local.py::test_complex_engine_does_not_call_global_redraw_candidates` | Clean2D | Engine reports ok where test expected preserved/no-op complex result. |

## Files eliminated or moved

| File eliminated/moved | Motive | Evidence | Replacement/coverage |
| --- | --- | --- | --- |
| `src/main.py` deleted | Compatibility stub for old `python src/main.py` execution. | Official entry point is `chemuson = "chemuson.__main__:main"`; docs were the only live references. | README and manual now document `chemuson`. |
| `src/verify_rendering.py` deleted | Manual visual/debug script living in `src`, with print/assert workflow and temp file side effects. | No imports from src/tests/tools/packaging; overlapped existing rendering/persistence tests. | Existing render/export/template/persistence tests remain. |
| `src/debug_h_angles.py` deleted | Local debug script with top-level QApplication/canvas side effects. | No integration with package, docs, CI, or tests. | Hydrogen rendering behavior remains covered by UI/model tests. |
| `repro_rdkit.py` deleted | One-off RDKit reproduction script. | No imports or docs/CI references; RDKit roundtrip/validation suites cover isolated worker and conversion behavior. | `tests/test_rdkit_roundtrip.py`, `tests/test_rdkit_smoke.py`, `tests/test_rdkit_strict_validation.py`. |
| `src/verify_orbitals_rendering.py` moved to `tools/visual/verify_orbitals_rendering.py` | Useful opt-in visual regression helper, but not product package code. | No package importers; has CLI/baseline comparison value. | Tool preserved under `tools/visual/`; path root fixed for repo layout. |
| `src/chemuson/gui/canvas/_shared.py` deleted | Legacy canvas facade with no real consumers. | Grep found no imports of `chemuson.gui.canvas._shared` or local canvas `_shared`; similarly named command `_shared` is separate and still live. | Canvas modules import real constants/data directly. |
| `src/chemuson/gui/main_window.py` updater reexports removed | Legacy API surface only used by one test. | Grep found `FLATPAK_APP_ID`, `format_no_update_message`, `format_update_disabled_message` imported from `main_window` only in `tests/test_update_ui_text.py`. | Test now imports from `chemuson.gui.controllers.update_controller`. |

## Tests eliminated or consolidated

| Tests eliminated/consolidated | Motive | New coverage location |
| --- | --- | --- |
| `tests/test_chemname_pr0.py` | PR-numbered historical wrapper around core `iupac_name(None)` behavior. | `tests/test_chemname_core.py::test_iupac_name_returns_not_available_for_empty_input`. |
| `tests/test_chemname_pr1.py` | PR-numbered MolView acyclic/cyclic checks with duplicated helper boilerplate. | `tests/test_chemname_core.py::test_molview_detects_acyclic_chain` and `test_molview_detects_cycle`. |
| `tests/test_template_match_prA.py` | PR-lettered template matching tests. | `tests/test_chemname_special_templates.py`. |
| `tests/test_chemcalc_pr6.py` | PR-numbered formula/mass tests. | `tests/test_chemcalc_formula.py`. |
| `sys.path.append/insert` blocks across tests | Redundant with pytest configuration and obscured real import roots. | Removed from all tests; `pyproject.toml` now uses `pythonpath = [".", "src"]` so tests can import both package code and `tools`. |
| Unused imports in tests after path cleanup | Dead import noise from removed `os`/`sys` blocks and historical imports. | `ruff check tests --select F401 --fix` removed 250 unused imports. |

## Functions/classes/imports eliminated

| Item | Motive | Evidence |
| --- | --- | --- |
| 28 unused imports in src/packaging/tests | Dead import cleanup. | `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` reported them; final run passes. |
| `Bond` annotation fix in `rdkit_io.py` | Existing lint failure after checking selected rules. | `F821` for `Bond` annotations; imported `Bond` from model. |
| Ambiguous `l` variable in `clean2d/safety.py` | Existing `E741` lint failure. | Renamed generator variable to `length`; no behavior change. |

## Commands executed

| Command | Final result |
| --- | --- |
| `python -m compileall src tests tools packaging` | Passed. |
| `pytest -q` | 842 passed, 55 skipped, 6 failed in ~50s. Same Clean2D failure set as baseline. |
| `pytest --collect-only -q` | 903 tests collected in ~0.8s. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Passed. |

## Remaining risks

| Risk | Notes |
| --- | --- |
| Clean2D baseline is red | Six Clean2D tests failed before cleanup and still fail. This cleanup intentionally did not change Clean2D behavior. |
| Many ChemName PR-numbered files remain | Only the clearly trivial/semantic PR wrappers were consolidated in this pass. Remaining PR files still encode chemistry cases and need domain-aware consolidation, likely into acceptance YAML or semantic suites. |
| `tools/visual/verify_orbitals_rendering.py` is opt-in | It is preserved but not wired into CI; visual regression remains manual unless a future opt-in marker/job is added. |

## Recommended next phase

1. Fix or explicitly rebaseline the 6 Clean2D failures before deeper Clean2D pruning.
2. Continue ChemName consolidation by moving repetitive PR cases into `tests/data/chemname_acceptance_cases.yml` or semantic suites.
3. Add a lightweight CI/pytest marker strategy for opt-in visual checks under `tools/visual` if orbital visual regression should be automated.
