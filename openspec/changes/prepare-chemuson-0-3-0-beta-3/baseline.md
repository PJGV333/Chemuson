# Baseline previo a los cambios de producción

- Fecha UTC: 2026-10-25.
- Rama: `release/v0.3.0-beta.3-prep` (creada desde `6aecba5a5585d9560b693186578b76bea9e623b8`, HEAD de `release/v0.3.0-beta.2-prep` y commit objetivo del tag `v0.3.0-beta.2`).
- `git status --short` al capturar: vacío (sin cambios; el scaffold OpenSpec de esta tarea se creó después de la captura de `baseline-git-status.txt`, que quedó vacío).

Las salidas en los archivos `baseline-*.txt` adjuntos son copia exacta de stdout/stderr capturado para cada comando. `baseline-pytest-collect.txt` conserva el listado completo de los 2639 tests.

## Comandos y códigos de salida

| Comando | Exit | Resultado |
|---|---:|---|
| `git status --short` | 0 | Sin cambios. Salida exacta: `baseline-git-status.txt` (vacío). |
| `python -m compileall src tests tools packaging` | 0 | Correcto. Salida exacta: `baseline-compileall.txt`. |
| `pytest --collect-only -q` | 0 | 2639 tests collected in 0.99s. Salida completa exacta: `baseline-pytest-collect.txt`. |
| `pytest -q` | No ejecutado | Omitido por la restricción explícita de no ejecutar la suite completa. Anotación exacta: `baseline-pytest-full.txt`. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | 1 | Misma única excepción basal heredada del baseline beta.2: F401 `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. No se corrige por estar fuera de alcance. Salida exacta: `baseline-ruff.txt`. |

Nota: la colección base era 2630 tests en beta.2; 2639 incluye ya las 9 pruebas nuevas de regresión del preflight añadidas por este cambio (la captura de colección se realizó tras añadirlas, antes de tocar workflows/versión; `baseline-git-status.txt`, `baseline-compileall.txt` y `baseline-ruff.txt` se capturaron sobre el árbol inmaculado de la rama).
