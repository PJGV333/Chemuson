# Baseline — add-editable-3d-molecular-editor

Fecha: 2026-10-16. Capturado antes de crear este OpenSpec o modificar código.

- Rama: `feature/chemuson-3d`, igual a `origin/feature/chemuson-3d`; `HEAD` `5ace94cf5ba97a14624452091c4840842a17d0a4`; `git status --short` vacío (exit 0).
- Entorno del proyecto: `.venv`, CPython 3.13.14, pytest 9.1.1, Ruff 0.17.1; RDKit local 2026.09.1.
- Se ejecutaron secuencialmente los checks del protocolo AGENTS. Por la discrepancia del Python del sistema, se usó el intérprete ejecutable del `.venv` del repositorio.

| Comando | Resultado |
|---|---|
| `git status --short` | exit 0, salida vacía |
| `./.venv/bin/python -m compileall src tests tools packaging` | exit 0 |
| `./.venv/bin/pytest --collect-only -q` | exit 0; 2102 tests collected in 0.45s |
| `./.venv/bin/pytest -q` | exit -11; avanzó aproximadamente al 61%, mostró dos marcadores `F` y terminó `Fatal Python error: Segmentation fault` en `concurrent.futures.thread._worker` con `<invalid frame>`. El log no contiene resumen/nombres porque el proceso murió. Fallo baseline conocido; no se investigará globalmente en esta campaña ni se atribuye al editor. |
| `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1; `F401` preexistente: `math` sin uso en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. |

Las salidas exactas (stdout/stderr capturados) están en `baseline/01-git-status.log` a `baseline/05-ruff.log`; `baseline/04-pytest.log` preserva el traceback del SIGSEGV. No se alteró ni se filtró el baseline. Los archivos de salida de colección suman aproximadamente 209 kB porque conservan la lista exacta de tests.
