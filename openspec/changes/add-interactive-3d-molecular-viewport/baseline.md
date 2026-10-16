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
- `timeout 600 pytest -q`: **exit -11** tras 281s; se completaron 1296 resultados (**1285 passed, 9 skipped, 2 failed**) antes de `Fatal Python error: Segmentation fault` en worker `concurrent.futures`/frame principal inválido.
- Los fallos completos fueron `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs` (#847) y `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` (#1074). No se atribuyen a este cambio.
- Ruff `src tests tools packaging --select F401,F811,F821,E722,E741`: **exit 1**, un `F401` preexistente para `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Salida completa de esta captura guardada fuera del repositorio en `/home/ccachyavgp/.hermes/cache/scratch/chemuson-viewport-stabilization-baseline-20261015.txt`.

## Comparación tras la estabilización

- `pytest --collect-only -q`: **2102 tests collected**, exit 0.
- `pytest -q`: **exit 139**; se completaron 1362 resultados (**1353 passed, 7 skipped, 2 failed**) antes del mismo patrón de segfault en worker `concurrent.futures`/frame principal inválido. El proceso se interrumpió antes de ejecutar las pruebas `viewport3d`.
- Los dos `F` aparecen en los mismos IDs e índices que en la captura pre-cambio (#847 y #1074). El primero también falla ejecutado aisladamente, por una expectativa de origen `rdkit_isolated`; el segundo pasa aisladamente y sólo falla en la suite completa. Se dejaron sin cambios y sin atribución causal.
- Las pruebas focalizadas de viewport, geometría 3D, CompChem y arquitectura: **154 passed, 1 skipped**. Compilación, OpenSpec estricto y Ruff del alcance: correctos. Ruff global conserva únicamente el `F401` ya documentado.
- La prueba de interacción real reportó GL listo, escenas de etano/benceno de **8/7** y **12/12** átomos/enlaces, picking del átomo visible `1`, orbit/pan/zoom efectivos y limpieza de la escena efímera al cambiar de documento.
- El proceso de interacción y captura guardó imágenes, pero terminó con `SIGSEGV` durante el apagado Qt/SIP. El mismo fallo se reprodujo al cerrar un `ChemusonWindow` con dos documentos y el dock 3D oculto, sin cargar escena ni inicializar el renderer. Esto demuestra que el fallo de apagado no requiere la ruta del visor, pero no determina la causa raíz; queda sin atribución.
- Capturas: `/home/ccachyavgp/.hermes/cache/scratch/chemuson_viewport_ethane_stabilized.png` y `/home/ccachyavgp/.hermes/cache/scratch/chemuson_viewport_benzene_stabilized.png`. Se verificó visualmente el anillo bencénico, los hidrógenos y las varillas; las capturas son evidencia del render, no de un apagado limpio.
