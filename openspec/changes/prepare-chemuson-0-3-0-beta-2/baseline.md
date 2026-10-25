# Baseline previo a los cambios de producción

- Fecha UTC: 2026-10-25.
- Rama: `release/v0.3.0-beta.2-prep` (aislada desde `release/v0.3.0-beta.1-prep`).
- SHA de código: `aa2e91c98cf1a16e8299e881955b1614edb2e18d`.
- El único estado no rastreado al capturar `git status --short` era el scaffold OpenSpec de esta tarea; no había cambios de producción.

Las salidas en los archivos `baseline-*.txt` adjuntos son copia exacta de stdout/stderr capturado para cada comando. `baseline-pytest-collect.txt` conserva el listado completo de los 2630 tests (no sólo el resumen).

## Comandos y códigos de salida

| Comando | Exit | Resultado |
|---|---:|---|
| `git status --short` | 0 | Sin cambios de producción; sólo el directorio OpenSpec de esta tarea estaba sin rastrear. Salida exacta: `baseline-git-status.txt`. |
| `python -m compileall src tests tools packaging` | 0 | Correcto. Salida exacta: `baseline-compileall.txt`. |
| `pytest --collect-only -q` | 0 | 2630 tests collected in 0.63s. Salida completa exacta: `baseline-pytest-collect.txt`. |
| `pytest -q` | No ejecutado | Omitido por la restricción explícita de no ejecutar la suite completa. La anotación exacta está en `baseline-pytest-full.txt`. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | 1 | Una excepción basal F401: `math` importado pero no usado en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. No se corrige por estar fuera de alcance. Salida exacta: `baseline-ruff.txt`. |

La identidad completa y salida de estado se conservan en `baseline-git-status.txt`.
