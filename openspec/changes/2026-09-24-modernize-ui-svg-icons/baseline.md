# Baseline — 2026-09-24-modernize-ui-svg-icons

Capturada el 2026-09-24 en la rama `ui/modernization` (HEAD `f1264fd`), **antes**
de modificar código, con el entorno reproducible de la Fase 1
(`uv run --no-project --offline` sobre el intérprete del venv principal).

## Entorno de validación

```
uv run --no-project --offline \
  --python /home/unison-pjgv/Documentos/GitHub/Chemuson/.venv/bin/python \
  --with pytest --with ruff --with PyQt6 --with numpy --with Pillow \
  --with rdkit --with certifi --with PyYAML -- <cmd>
```

## git status --short

(árbol limpio, sin salida)

```


```

## python -m compileall src tests tools packaging

```
compileall rc=0

```

## pytest --collect-only -q (resumen)

```
tests/test_wedge_geometry.py::test_aromatic_double_bold_uses_thin_secondary_pi_line

1636 tests collected in 0.61s
```

## pytest -q (suite completa)

```
=========================== short test summary info ============================
FAILED tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs
FAILED tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash
FAILED tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo
4 failed, 1612 passed, 20 skipped in 330.73s (0:05:30)
pytest rc=1
```

**Lectura**: 1612 passed / 20 skipped / 4 failed. Los 4 fallos son
**preexistentes e idénticos** a la baseline de la Fase 1
(`2026-09-24-modernize-ui-theme-foundation/baseline.md`):
- `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`
- `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`
- `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`
- `tests/test_smiles_stereo_import.py::test_tetrandrine_import_preserves_visual_stereo`

Se dejan sin tocar (fuera de alcance; AGENTS.md §5).

## ruff check src tests tools packaging --select F401,F811,F821,E722,E741

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
  |
2 |
  - import math
3 |
  |

Found 1 error.
[*] 1 fixable with the `--fix` option.
ruff rc=1

```

**Lectura**: 1 finding **preexistente** (F401 `import math` en
`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`), idéntico al
de la Fase 1; fuera de alcance, no se corrige.

## Conclusión

Baseline estable y equivalente a la de la Fase 1. Criterio de parada de
AGENTS.md: no se aplica (no hay fallos nuevos inesperados ni dependencias
no catalogadas).
