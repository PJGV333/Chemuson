# Baseline — antes de cualquier modificación de código

Rama: `fix/primary-amine-label-orientation-beta1` creada desde
`origin/stabilize/v0.3.0-beta.1` @ `33ed522381dce83fcc011e4fb0241919440ef32f`.

## `git status --short`

Salida vacía (worktree limpio antes de modificar).

## `chem/bin/python -m compileall -q src tests tools packaging`

PASS (sin salida de errores).

## `pytest -q` (suite completa, entorno local Python 3.14.8 / rdkit 2026.03.6)

**63 failed, 2547 passed, 22 skipped** en 9:45.

Clústeres de fallo (preexistentes en HEAD antes de tocar el código):
`tests/test_aromatic_nitrogen_regression.py` (numerosos),
`tests/test_chemname_smoke.py::test_rdkit_descriptors_smoke`,
`tests/test_functional_group_hydrogen_adjustment.py` (varios),
`tests/test_implicit_h_substitution.py` (3). En la campaña anterior la suite verde en
CI cubrió los 2.610 node IDs; estos fallos no aparecen al ejecutar esos mismos archivos
de forma aislada (antes y después de este cambio), por lo que se clasifican como
ruido de entorno/interferencia entre tests de la corrida monolítica local
(ver `validation.md`). No se modificó ninguno de esos tests como parte de este cambio.

## `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

Un único hallazgo preexistente: F401 `math` sin usar en
`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py` (documentado en
campañas previas; fuera de alcance, no se tocó).
