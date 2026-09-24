===== git status --short =====
?? openspec/changes/2026-09-24-modernize-ui-theme-foundation/
===== EXIT 0 =====

===== git status --short --branch =====
## ui/modernization...origin/ui/modernization
?? openspec/changes/2026-09-24-modernize-ui-theme-foundation/

===== python -m compileall src tests tools packaging =====
EXIT 0

===== pytest --collect-only -q =====
ERROR tests/architecture/test_main_window_clean2d_geometry.py
ERROR tests/architecture/test_module_catalog.py
ERROR tests/architecture/test_platform_settings_boundary.py
ERROR tests/architecture/test_public_api_exists.py
ERROR tests/architecture/test_recovery_boundary.py
ERROR tests/architecture/test_update_boundary_audit.py
!!!!!!!!!!!!!!!!!!! Interrupted: 16 errors during collection !!!!!!!!!!!!!!!!!!!
1399 tests collected, 16 errors in 5.43s


===== pytest -q (suite completa) =====
=========================== short test summary info ============================
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
4 failed, 1583 passed, 20 skipped in 399.66s (0:06:39)

===== Fallos preexistentes (baseline) =====
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo

===== ruff check src tests tools packaging --select F401,F811,F821,E722,E741 =====
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
  |
1 | from __future__ import annotations
2 |
3 | import math
  |        ^^^^
4 |
5 | from chemuson.clean2d import (
  |
help: Remove unused import: `math`
  |
2 |
  - import math
3 |
  |

Found 1 error.
[*] 1 fixable with the `--fix` option.

===== Notas del entorno de test =====
Worktree sin .venv. Se usa entorno efimero:
  uv run --no-project --offline --python /home/unison-pjgv/Documentos/GitHub/Chemuson/.venv/bin/python
  --with pytest --with ruff --with PyQt6 --with numpy --with Pillow --with rdkit --with certifi --with PyYAML
Intérprete: Python 3.11.16 (venv del checkout principal, stack runtime PyQt6 6.11.0/Qt 6.11.2).
pytest 9.1.1 y ruff resueltos desde la caché de uv (offline; sin instalaciones nuevas).
Fallos/ruff preexistentes: 4 tests (candidate generation y stereo import) + 1 F401 en
tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py (campaña Clean2D aislada, fuera de alcance).

=====================================================================
VERIFICACION FINAL (post-cambio, mismo comando que la baseline)
=====================================================================

===== git status --short (post) =====
 M AGENT_REPORT.md
 M architecture/modules.yml
 M src/chemuson/gui/main_window.py
 M src/chemuson/gui/shell/assembly.py
 M src/chemuson/gui/styles.py
 M src/chemuson/platform/__init__.py
 M src/chemuson/platform/settings.py
 M tests/test_platform_settings.py
?? docs/ui-modernization/foundation-shots/
?? openspec/changes/2026-09-24-modernize-ui-theme-foundation/
?? src/chemuson/gui/theme/
?? tests/test_ui_theme_foundation.py

===== python -m compileall src tests tools packaging =====
EXIT 0 (sin salida = sin errores)

===== pytest -q (suite completa, post) =====
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
4 failed, 1612 passed, 20 skipped in 335.48s (0:05:35)

===== Comparacion de fallos contra baseline =====
Identicos: SI - solo fallos preexistentes
Delta tests: +29 nuevos (28 en test_ui_theme_foundation.py + 1 en test_platform_settings.py)

===== ruff check (post) =====
1 error: el mismo F401 preexistente de tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py (sin nuevos errores)

===== git diff --check =====
OK (sin errores)

===== Smoke Qt offscreen (ventana real, light+dark) =====
SMOKE: OK — capturas en docs/ui-modernization/foundation-shots/
