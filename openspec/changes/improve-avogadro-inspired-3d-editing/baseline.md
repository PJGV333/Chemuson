# Baseline — improve-avogadro-inspired-3d-editing

Captura local: 2026-10-16 (zona local del entorno, UTC−07:00), antes de cambios de producto/documentación. Worktree: `/home/ccachyavgp/Documentos/Chemuson-3D`, independiente del worktree `/home/ccachyavgp/Documentos/ChemUSON`.

- Rama: `feature/chemuson-3d`; upstream: `origin/feature/chemuson-3d`.
- HEAD: `bfccb68c5e60221369f1ba83bfe701a4df3e8f8e`.
- `git status --short --branch`: `## feature/chemuson-3d...origin/feature/chemuson-3d` y `?? PROMPT_LUNA_TEMP.txt`. Ese archivo contiene el encargo temporal, no forma parte del producto y no se incluirá en ningún commit.
- `git worktree list --porcelain` confirmó que el worktree objetivo es la ruta Chemuson-3D en `feature/chemuson-3d`; el otro worktree está en otro repo/rama y queda intacto.
- `./.venv/bin/python -m compileall -q src tests tools packaging`: exit 0, salida vacía.
- `./.venv/bin/pytest --collect-only -q`: **2182 tests collected in 0.45s**, exit 0. Salida completa: `/home/ccachyavgp/.hermes/cache/scratch/chemuson-campaign-20261017-collect.log`.
- `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por el único diagnóstico preexistente `F401` (`math` sin uso) en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- `./.venv/bin/pytest -q`: **no ejecutado** en esta captura. Se conserva el baseline anterior de la rama: SIGSEGV/exit 139 en worker `concurrent.futures` durante la suite monolítica; consultar `openspec/changes/stabilize-3d-phase-1-4-integrity/validation.md` y su `baseline.md`. El usuario pide evitar reruns innecesarios; la omisión no cuenta como pase.

No se ejecutaron GUI, binarios de Avogadro, comandos privilegiados, ni cambios de sistema. No se hizo push, merge ni PR.
