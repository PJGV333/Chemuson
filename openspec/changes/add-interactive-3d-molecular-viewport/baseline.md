# Baseline — add-interactive-3d-molecular-viewport

Fecha de captura: 2026-10-15.

- Rama/estado: `feature/chemuson-3d...origin/feature/chemuson-3d`; árbol limpio antes de crear este cambio.
- HEAD observado: `f0fde60 chore(release): prepare v0.3.0-beta.1 and preview builds`.
- Entorno de pruebas: `.venv` creado por `uv`, CPython 3.13.14, PyQt6 6.11.0, RDKit disponible; `PyYAML`, pytest y Ruff instalados sólo en ese entorno de validación local. No se modificaron dependencias runtime del proyecto.
- `python -m compileall src tests tools packaging`: **exit 0**.
- `pytest --collect-only -q`: **2072 tests collected in 0.46s**, **exit 0**.
- `pytest -q`: **exit 139**. La ejecución alcanzó aproximadamente el 62 %, mostró fallos `F` preexistentes y terminó con `Fatal Python error: Segmentation fault` en un worker `concurrent.futures`. No se atribuye este fallo al cambio porque ocurrió antes de implementar el visor.
- Ruff `src tests tools packaging --select F401,F811,F821,E722,E741`: **exit 1** por el hallazgo preexistente `F401` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`).
- El shell no interactivo emitió advertencias `tcsetattr`/control de trabajos al finalizar pytest; se consideran limitación del entorno, no evidencia de éxito de la suite.

Salida completa conservada fuera del árbol de producción en:
`/home/ccachyavgp/.hermes/cache/scratch/chemuson-3d-viewport-baseline-20261015-complete.txt`

## Baseline de estabilización (antes de cambios, HEAD `e2a5a3d`)

- Fecha: 2026-10-15. Rama `feature/chemuson-3d`, `HEAD` igual a `origin/feature/chemuson-3d`; worktree limpio.
- `python -m compileall src tests tools packaging`: **exit 0**.
- `pytest --collect-only -q`: **2084 tests collected in 0.43s**, **exit 0**.
- `timeout 600 pytest -q`: **exit -11** tras 281s, progreso hasta aproximadamente el 62 %, dos `F` visibles al 41 % y 51 %, seguido por `Fatal Python error: Segmentation fault` con worker `concurrent.futures` y frame principal inválido. No se conservan nombres ni trazas de esos dos fallos; ninguno se atribuye a este cambio sin reproducción diferencial.
- Ruff `src tests tools packaging --select F401,F811,F821,E722,E741`: **exit 1**, un `F401` preexistente para `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Salida completa de esta captura guardada fuera del repositorio en `/home/ccachyavgp/.hermes/cache/scratch/chemuson-viewport-stabilization-baseline-20261015.txt`.
