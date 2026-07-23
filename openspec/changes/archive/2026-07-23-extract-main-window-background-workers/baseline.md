# Baseline

Fecha: 2026-07-23

- Rama inicial: `main` en `05bfa360`.
- `git status --short`: limpio.
- `python -m compileall src tests tools packaging`: correcto.
- `pytest --collect-only -q`: no ejecutable; `pytest` no está instalado en la
  imagen actual.
- `pytest -q`: no ejecutable por la misma limitación.
- Ruff: no ejecutable; `ruff` no está instalado en la imagen actual.
- `openspec validate --all --strict`: 20/20 válidos.

La fase anterior registró 210 pruebas focalizadas/arquitectónicas correctas y
el usuario confirmó manualmente el arranque y funcionamiento de Chemuson en
Qt real. La suite completa de aquella fase conservó 56 bloqueos ambientales
por ausencia de `libEGL.so.1`.
