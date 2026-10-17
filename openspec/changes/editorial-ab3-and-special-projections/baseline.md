# Baseline — Preset editorial AB3 y proyecciones especiales

## Estado del repositorio

- Rama: `feat/molecular-depiction-visual-lab`.
- HEAD: `70efde4bca68adc6ffc55a10bd18727fb6fa1cec`.
- `git status --short`: vacío antes de iniciar la propuesta; HEAD limpio es el checkpoint previo a cambios.
- Transcripción completa de baseline: `baseline-command-output.txt`.

## Comandos obligatorios

| Comando | Resultado basal |
|---|---|
| `python -m compileall src tests tools packaging` | PASS, exit 0. |
| `pytest --collect-only -q` | PASS, 2,326 tests collected in 0.72 s. |
| `pytest -q` | 2,305 passed, 20 skipped, 1 failed in 1,256.22 s. El fallo es `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`, por orden/lifecycle Qt (`worker.finished` observado sin `thread.finished`/`job_finished`). Es el mismo fallo heredado y dependiente del orden ya documentado en el OpenSpec previo del laboratorio; se conserva como baseline, no se atribuye a este cambio. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el F401 preexistente `math` sin uso en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; fuera de alcance. |

El SHA y la transcripción se registraron antes de editar código o documentación. La ejecución basal monolítica excedió 20 minutos; el cierre final usará las particiones acotadas del CI, reconciliadas contra colección total, y reportará el fallo heredado sin ocultarlo.
