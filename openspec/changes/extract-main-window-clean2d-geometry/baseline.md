# Baseline

Fecha: 2026-07-23

- Rama inicial: `architecture/phase9-extract-main-window-clean2d-geometry`,
  basada en `main` `97a3f032`.
- `git status --short`: limpio antes de archivar la fase 8.
- `python -m compileall src tests tools packaging`: correcto.
- El Python predeterminado no incluye pytest ni Ruff. El entorno persistente
  `/workspace/.venvs/chemuson` sí los incluye.
- `pytest --collect-only -q` desde ese entorno: 1050 pruebas recolectadas y 56
  errores de importación por ausencia de `libEGL.so.1`.
- `pytest -q`: reproduce los mismos 56 errores durante colección.
- Ruff global: reproduce únicamente el F401 histórico de
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.
- `openspec validate --all --strict`: 21/21 válidos antes y después de
  archivar la fase 8.

La fase 8 fue probada manualmente en Qt real y aprobada por el usuario antes
de su archivado.
