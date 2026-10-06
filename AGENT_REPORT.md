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

## Continuación ChemUSON — Phase 4.5/5 (2026-10-07)

- Rama `ai/molecular-assistant-foundation`, base `6d4fca961cc7695ed01d78ae8f0a1ad769388384`. `c0233eb` (funcional) y `6d36e3a` (higiene/dependencias) se publicaron con pushes normales. No hubo merge, rebase ni force-push. El commit de cierre de este informe/tasks se publicará también mediante push normal.
- Se cerró el contrato tipado de transformación M23, Insert Variant/Replace undoable, controles y perfiles runtime sin secretos, verificación conservadora de identidad con ChemIO aislado y resumen seguro del evaluador Clean2D. No se modificó `src/chemuson/clean2d/`, persistencia `.cmsn` ni dependencias runtime.
- Tests: colección actual 1946; shards cubrieron 1922 passed, 20 skipped, 4 fallidos. Los fallos son el test Clean2D ya registrado en baseline y tres aserciones de importación estereoquímica bajo RDKit 2026.03.6. Suite de arquitectura 278 passed; OpenSpec estricto 52 passed. Ruff de archivos cambiados PASS; Ruff global conserva solo el F401 `math` histórico.
- No se repitió la suite monolítica: su referencia histórica es 19:26 y el límite activo es 10 minutos. El par CompChem→Assistant pasó 5/5 en ambos órdenes; el bloque ordenado Assistant→transform→UI pasó 99. No se reprodujo SIGSEGV en los shards, pero tampoco se afirma que el aborto monolítico histórico esté descartado ni se encontró un segundo trigger de producción.
- Smoke local únicamente contra Qwen 3.8 27B en `127.0.0.1:8081`: etanol/cafeína pasaron ChemIO e identidad offline; colesterol, vancomicina y eritromicina devolvieron `invalid_json`. No se usó red externa ni se guardó salida cruda. Evidencia detallada: `openspec/changes/stabilize-ai-molecular-assistant-integration/validation.md` y `/tmp/chemuson-qwen-local-smoke/`.
