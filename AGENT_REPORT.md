# AGENT_REPORT — cierre de higiene del repositorio

## Alcance y decisiones

- Rama `maintenance/repository-hygiene-closure`, creada desde `origin/main`
  `1db4f63b52af79247745b3a8a220fb728348218c`, publicada mediante push normal.
- `origin/main` permanece en el mismo SHA. No se hizo merge, force-push,
  rebase, cherry-pick, `git gc --prune` ni reescritura de historia.
- OpenSpec: `openspec/changes/repository-hygiene-closure-2026-10-03/`;
  checklist completo. Memoria y decisiones en
  [`docs/history/CAMPAIGNS.md`](docs/history/CAMPAIGNS.md), política en
  [`docs/history/REPOSITORY_POLICY.md`](docs/history/REPOSITORY_POLICY.md) y
  métricas/refs/prune en
  [`docs/history/REPOSITORY_CLEANUP_2026-10-03.md`](docs/history/REPOSITORY_CLEANUP_2026-10-03.md).

## Cambios

- Se retiraron `src/sys`, `src/repro_v2.png`, dos parches sueltos ya históricos,
  124 outputs orbitales bajo `tests/archive`, prototipo UI ejecutable y
  capturas/generadores intermedios, scripts F7 one-shot y cuatro informes
  redundantes. OpenSpec Markdown/specs, tests y fixtures vigentes, cinco
  capturas UI finales, tres capturas KDE/Wayland, siete evidencias F7 selectas
  y el `theme.py` normativo se conservan.
- Se corrigieron enlaces absolutos `file://` del manual, se acortaron las
  fichas UI/known issues y comentarios de tema que apuntaban a prototipos
  retirados. Ningún archivo test de código ni módulo `src/chemuson/clean2d/`
  fue modificado; tampoco se cambiaron química, geometría, serialización o
  versión.
- Los previews orbitales ahora escriben fuera del árbol por defecto. El
  `render_orbital_family_preview.py` pasó smoke; el `orbital_fit_report.py`
  falla por un bug ya existente (`PiBondingParams.ring` inexistente en
  `_family_metric_strings`), dejando salidas parciales solo en temp. Se
  documenta y no se arregla fuera de alcance.
- Se borraron 23 ramas remotas y tres locales solo tras verificar SHA y
  ancestro de `origin/main`; quedaron intactas ocho ramas únicas/activas,
  incluida `gh-pages` y ambas Clean2D. Inventario final: 10 remotas reales,
  cuatro locales.

## Verificación

- `python -m compileall src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 1816.
- `pytest -q tests/architecture`: 269 passed.
- UI dirigida: 299 passed.
- `pytest -q`: 1760 passed, 55 skipped, 1 fallo de baseline persistente:
  `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`.
  Fue el único fallo en baseline y en ambos full runs del cierre; no se cambió.
- Ruff scoped `F401,F811,F821,E722,E741`: un F401 histórico en el test
  Clean2D `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; no se
  toca fuera de alcance.
- `openspec validate --all --strict`: 43 passed, 0 failed.
- Smoke Qt offscreen 100%/HiDPI 200%: PASS (temas claro/oscuro, 1440×900,
  980×600, DPR 2, tabs, rail/flyout, SidePanel, paleta y onboarding).
- `git diff --check`: PASS; escaneo de 350 Markdown: cero links relativos
  rotos. `git status --short`: limpio después del commit/push final.

Una primera invocación del smoke Qt olvidó `PYTHONPATH=src`; se repitió con
`PYTHONPATH` y `XDG_CONFIG_HOME` temporales y pasó, sin tocar preferencias del
usuario. Un primer `diff --check` encontró tres espacios finales en el baseline
OpenSpec; se corrigieron antes del push.

## Integridad final

Los tamaños, el alcance exacto de archivos, SHAs previos/posteriores y la
estimación de compactación histórica están en el reporte enlazado arriba. Los
objetos Flatpak activos publicados por `gh-pages` se preservan. No se inició
otra fase de producto.
