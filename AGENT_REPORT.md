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

## Bloqueo de diseño — Phase 5 Molecular Assistant

**Estado al `ed2953d728c3dad78ad3adb085b89256bbb5d489`** (`ai/molecular-assistant-foundation`, limpio y sincronizado con `origin`). Phase 2, 3 y 4 están implementadas, validadas y publicadas en commits separados. No se hizo merge a `main`.

No se implementa Phase 5 todavía: “edición molecular estructurada mediante IA” no fija la semántica de edición ni la forma del cambio, y esas alternativas afectan operaciones químicas y comportamiento undo/redo incompatibles. El repositorio ofrece `molgraph_to_smiles_isolated`, validación M23 de SMILES completo, `DeleteSelectionCommand` y la inserción undoable del canvas; aun así hay que elegir entre:

1. **Reemplazo completo del componente molecular seleccionado**: permitir exactamente una molécula aislada completa; exportarla a SMILES aislado, enviar la instrucción y el SMILES fuente al M23 existente, mostrar fuente/propuesta y reemplazarla sólo tras aprobación explícita, dentro de un único paso undoable. La salida se conserva como SMILES completo y pasa por la validación ChemIO ya existente. No habría DSL ni cambios al contrato público de M23.
2. **Operaciones atómicas estructuradas**: devolver y aplicar una lista allowlisted de acciones sobre átomos/enlaces. Esto requiere decidir operaciones exactas, cómo se identifican átomos sin ambigüedad frente al mapeo SMILES/MolGraph, qué restricciones de valencia/estereoquímica aplicar y la política de rechazo/undo.
3. **Propuesta no destructiva como estructura separada**: insertar una molécula derivada nueva y conservar la fuente intacta; es más segura, pero no reemplaza/edita la selección.

La opción 1 parece el alcance mínimo que conserva el contrato M23 aprobado y entrega edición real con confirmación y undo; no se adopta sin autorización porque cambia la semántica de la selección existente. **Decisión solicitada:** confirmar la opción 1, elegir la opción 3 o definir el allowlist/mapeo de la opción 2. Hasta entonces no se crea el OpenSpec normativo de Phase 5 ni se modifica código/canvas/contratos químicos. No se introdujeron agentes autónomos, tool calling ni ejecución de código.
