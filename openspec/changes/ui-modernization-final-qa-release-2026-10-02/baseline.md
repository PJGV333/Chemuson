# Baseline — Fase 8 UI Modernization Final QA

Capturada en `release/ui-modernization-qa` antes de crear los artefactos OpenSpec de Fase 8 y antes de otros cambios.

## Estado y revisión

- `git status --short`: sin salida (worktree limpio).
- `git rev-parse HEAD`: `140a080c336515650bbeea0a4e6ead67a9999b23`.
- `git log -5 --oneline`:
  - `140a080 Normalize UI OpenSpec requirements after integration`
  - `546d08c Fix onboarding rendering and template click UX`
  - `5b1dc1d Fix template thumbnail HiDPI scaling`
  - `0dd4e0a Fix Fase 7 polish edge cases`
  - `3d4e4a8 Mark Fase 7 task 10.4 complete (commits + push verified)`
- Python: `3.14.7`.

## Compilación y colección

- `python -m compileall src tests tools packaging`: exit 0.
  Log completo de la captura: `/tmp/f8_compileall.log`.
- `pytest --collect-only -q`: exit 0; `1816 tests collected in 0.75s`.
  Log completo: `/tmp/f8_collect.log`.

## Suite completa

Comando: `pytest -q > /tmp/f8_full_suite.log 2>&1`.

- Exit: 1.
- Resultado: `1 failed, 1760 passed, 55 skipped in 1137.93s (0:18:57)`.
- Único fallo:
  `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`.
- Identidad verificada contra `openspec/changes/2026-10-01-modernize-ui-polish/baseline.md`: es uno de los cinco fallos preexistentes documentados allí como fallos RDKit/async. No se modifica el test ni el controller.
- Salida solicitada (tail -n 100) revisada; log íntegro en `/tmp/f8_full_suite.log`.

No se hicieron cambios de código durante la captura de baseline.