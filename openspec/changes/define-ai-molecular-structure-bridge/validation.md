# Validación — define-ai-molecular-structure-bridge

Validación final del cambio documental en `ai/molecular-assistant-foundation`, basado en `origin/main` `4068319e9a1deee7dbf239872191a8f509c8b641`. No hay cambios de runtime ni dependencias.

## OpenSpec y arquitectura

```text
$ openspec validate define-ai-molecular-structure-bridge --strict
Change 'define-ai-molecular-structure-bridge' is valid

$ openspec validate --all --strict
Totals: 44 passed, 0 failed (44 items)

$ pytest tests/architecture -q
269 passed in 9.24s

$ git diff --check
(exit 0)
```

## Verificación general

```text
$ python -m compileall src tests tools packaging
(exit 0)

$ pytest --collect-only -q
1816 tests collected in 0.75s
(exit 0)

$ pytest -q
FAILED tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend
1 failed, 1760 passed, 55 skipped in 1114.01s (0:18:34)
(exit 1)
```

La identidad del único fallo y los conteos coinciden con el baseline anterior a los cambios; es el fallo histórico descrito en `baseline.md`/`docs/history/CAMPAIGNS.md`, no una regresión de este cambio documental.

```text
$ ruff check src tests tools packaging --select F401,F811,F821,E722,E741
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
Found 1 error.
(exit 1)
```

El F401 es idéntico al baseline, preexistente y fuera de alcance; no se alteró.

## Alcance comprobado

- Sólo se crean artefactos OpenSpec de esta campaña.
- No se modifican `architecture/modules.yml`, código fuente, tests existentes, Clean2D, dependencias ni configuración del producto.
- La rama es independiente y parte del `origin/main` vigente; el push normal y el commit se registran en el informe de cierre de la sesión una vez realizados.