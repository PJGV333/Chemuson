# Baseline — integración 3D en visual lab Beta 1

## Identidad y estado inicial

| Elemento | Valor |
|---|---|
| Worktree | `/home/ccachyavgp/Documentos/Chemuson-3D-beta1-integration` |
| Rama local | `integration/chemuson-3d-beta1` |
| HEAD/base visual | `43a8a27c45cdabfd157bf3f2a52c9216b24a1552` (`docs(openspec): record remote visual lab validation`) |
| Rama fuente 3D | `origin/feature/chemuson-3d` = `ec486626333ccd2e0ae86671861708b30c508c5b` |
| Checkpoint local | `checkpoint/chemuson3d-beta1-base-43a8a27` en `43a8a27c45cdabfd157bf3f2a52c9216b24a1552` |
| Release observada, no modificada | `origin/release/v0.3.0-beta.1-prep` = `9a70ee1c908a39e700098bcc0a54d464a606e24c` |
| `git status --short` | Sin salida; árbol limpio antes de esta especificación |
| Historial | `git rev-list --left-right --count origin/feat/molecular-depiction-visual-lab...origin/feature/chemuson-3d` → `85 36` commits exclusivos respecto al ancestro común |

## Comandos y resultados exactos relevantes

- `python -m compileall -q src tests tools packaging` con el `python` por defecto del entorno Hermes: exit 0.
- `python -m pytest --collect-only -q` y las pruebas focales con ese `python` fallaron antes de iniciar porque ese intérprete no tiene instalado `pytest`. El worktree temporal no tiene `.venv` propio.
- Intérprete válido usado para pytest/compileall: `/home/ccachyavgp/Documentos/Chemuson-3D/.venv/bin/python`, Python **3.13.14**, pytest **9.1.1**. `pyproject.toml` añade `.` y `src` al `pythonpath`, así que pytest resuelve el código desde el worktree temporal.
- `/home/ccachyavgp/Documentos/Chemuson-3D/.venv/bin/python -m compileall -q src tests tools packaging` → exit 0.
- `/home/ccachyavgp/Documentos/Chemuson-3D/.venv/bin/python -m pytest --collect-only -q` → **2358 tests collected in 0.53s**.
- `/home/ccachyavgp/Documentos/Chemuson-3D/.venv/bin/python -m pytest -q tests/test_clean2d*.py tests/test_*depict*.py tests/architecture/test_main_window_clean2d_geometry.py` → **450 passed, 3 skipped in 170.94s**.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` → exit 1 por el hallazgo baseline `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificará fuera de alcance.
- `pytest -q` monolítico no ejecutado por la restricción explícita vigente de no repetirlo debido al SIGSEGV histórico Qt/SIP. Se evaluarán las regresiones críticas y el workflow particionado sin ocultar fallos ni timeouts.
- Previsualización read-only `git merge-tree --write-tree --name-only origin/feat/molecular-depiction-visual-lab origin/feature/chemuson-3d` → exit 1 esperado por conflictos; no modificó el índice ni el worktree. Conflictos: `AGENT_REPORT.md`, `architecture/modules.yml`, `src/chemuson/chemio/_rdkit_worker.py`, `src/chemuson/chemio/rdkit_io.py`, `src/chemuson/chemio/rdkit_safe.py`, `src/chemuson/gui/main_window_ui_builder.py`. Avisos de auto-fusión: `canvas_structure.py`, `canvas_view.py`, `docks.py`, `main_window.py`, `shell/assembly.py`, `tests/test_compchem3d_dock.py`.

## Limitaciones ya documentadas en la rama visual

Los OpenSpecs de depiction y preparación de Beta 1 contienen revisión visual del propietario pendiente y validaciones de suite/Qt con timeout o fallo dependiente del orden. La campaña de visualización de proyecciones también deja aceptación manual pendiente. Estos estados son gates, no se reinterpretarán como aprobados ni se borrarán al fusionar. Referencias: `openspec/changes/harden-molecular-depiction-and-drawing/validation.md`, `openspec/changes/refine-wedge-double-bond-junctions/validation.md`, `openspec/changes/prepare-chemuson-0-3-0-beta-1/validation.md` y `openspec/changes/stabilize-visual-lab-ci-and-projections/validation.md`.