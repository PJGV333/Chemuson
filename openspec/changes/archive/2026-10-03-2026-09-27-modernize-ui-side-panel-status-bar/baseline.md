# Baseline — 2026-09-27-modernize-ui-side-panel-status-bar

Capturada antes de modificar código o documentación de esta fase, en
`ui/modernization`, con árbol de trabajo limpio.

## Checkout

```
git status --short
(sin salida; árbol limpio)
git branch --show-current
ui/modernization
git rev-parse HEAD
d2468a3cd1b21aae70f1413c930c2a9bc46c599f
git rev-parse origin/ui/modernization
d2468a3cd1b21aae70f1413c930c2a9bc46c599f
```

## Entorno observado

- `python` / `python3`: CPython 3.11.15 (compilación).
- `/usr/bin/python`: CPython 3.14.7; pytest 9.1.1 y PyQt6 instalados (suite).
- Ruff 0.14.14.
- OpenSpec 1.5.0.
- Sin cambios de dependencias.

## `python -m compileall src tests tools packaging`

```
exit code: 0
```

## `pytest --collect-only -q`

```
1799 tests collected in 1.99s
```

## `pytest -q`

```
1744 passed, 55 skipped in 558.06s (0:09:18)
```

La suite de baseline no presentó fallos. Los cuatro fallos RDKit mencionados
como históricos en el encargo no se reprodujeron en este entorno; no se
modificó ningún test ni se hizo ningún arreglo incidental.

## `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

```
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

Found 1 error.
(exit code: 1)
```

Aviso preexistente, fuera del alcance de Fase 5; se deja intacto.

## Baseline previo al polish visual (2026-09-27)

Capturado en el árbol limpio de `ui/modernization`, HEAD local/remoto
`801e00376c4c28d3a9077508f9e6cff3ee0c067f`, antes de cambiar SideTabRow.

### Estado y verificaciones

- `git status --short`: sin salida; árbol limpio.
- `python -m compileall src tests tools packaging`: exit code 0.
- `pytest --collect-only -q`: `1811 tests collected in 0.76s`.
- `pytest -q` en este mismo HEAD: `1756 passed, 55 skipped in 756.44s (0:12:36)`.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`:

```
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

Found 1 error.
[*] 1 fixable with the `--fix` option.
```

El F401 es preexistente y no relacionado con el polish visual.
