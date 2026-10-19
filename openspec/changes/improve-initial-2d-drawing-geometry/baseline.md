# Baseline — geometría inicial 2D

## Revisión, rama y SHA

Se leyó completo el prompt de campaña en la raíz (se mantiene sin seguimiento y fuera de commits), el protocolo `AGENTS.md` y los OpenSpecs de depiction/junctions relacionados. `git fetch origin` confirmó que `origin/release/v0.3.0-beta.1-prep` contiene el merge `3e0ef606`; rama de trabajo `fix/initial-2d-geometry-beta1`, SHA base `3e0ef606edb90929c1fb32dddcd5385fabcc5d82`. No se trabajó en release, main, 3D o ramas Visual Lab.

### Línea de datos auditada

`canvas_selection_input.mousePressEvent` → `_begin_place_bond`/`_begin_place_ring`/`_begin_place_chain` → preview (`_update_*`) → helpers en `canvas_bond_geometry.py`/`canvas_tools_rings_chains.py` → release y comandos `AddBondCommand`/`AddRingCommand`/`AddChainCommand` → coordenadas en `MolGraph` → `AtomItem`/`BondItem.update_positions()` → `QGraphicsScene`/píxeles. `BondItem` crea paths recortados y offsets para dobles/aromáticos; no altera las coordenadas del modelo.

## Baseline de código y tests

Los stdout completos están guardados en `baseline/01`–`baseline/07`; el primero incluye la salida directa, y `baseline/06-first-failure.log.gz` se conserva comprimido sin alterar bytes para preservar una línea de pytest que trae whitespace final (`gzip -dc baseline/06-first-failure.log.gz`):

| Comando | Resultado basal |
|---|---|
| `git status --short` | Sólo `?? ChemUSON_Campana_Geometria_Inicial_2D_Beta1.md`; ningún cambio de usuario adicional |
| `python -m compileall src tests tools packaging` | PASS, exit 0 |
| `pytest --collect-only -q` | PASS, 2,567 tests, exit 0 |
| `timeout --signal=TERM 600s pytest -q` | Timeout a 600 s (exit 124), progresó hasta ~30%; mostró un `F` sin resumen por timeout |
| `timeout --signal=TERM 300s pytest -x -q` | Localiza el fallo basal: `tests/test_beta_compact_shell_and_tabs.py::test_view_menu_submenus_reuse_actions_and_remove_obsolete_auxiliary_option`, 366 passed/1 failed en 38.99 s. El test espera siete entradas y el merge beta actual incluye «Visor molecular 3D». No se toca por estar fuera de alcance. |
| Suites focales existentes: sp3, zigzag, canvas edges, junctions, wedge/double/labels, Haworth renderer | PASS, 56 passed en 6.94 s |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Un F401 basal no relacionado: `math` sin uso en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no corregir fuera de alcance |

La suite completa no se repitió a ciegas: el resultado basal ya excedió 10 minutos y `-x` identificó un fallo basal específico. La aceptación final usará suites focales y arquitectónicas con límites.

## Medidas iniciales reproducibles

`evidence/reproduce_initial_2d_baseline.py` construye benceno y Ph–CH₂–C(=O)–CH₂–CH₃ mediante gestos reales QTest de arrastre, sin Clean2D, y guarda JSON del modelo e imágenes producidas por el renderer real de escena. Se ejecutó con `PYTHONPATH=src QT_QPA_PLATFORM=offscreen`; el artefacto generado pertenece al baseline y no equivale a revisión manual en pantalla.

| Caso/medida antes | Previsto | Obtenido | Lectura |
|---|---:|---:|---|
| Benceno libre: seis lados | 40.0 | 39.999999999999986–40.000000000000036 | PASS; error máximo <4e-14 |
| Benceno libre: ángulo interior | 120° | 119.99999999999989–120.00000000000011° | PASS; error máximo <1.2e-13° |
| sp3 por sugerencia con un vecino (sin rejilla fija) | 120° editorial | 109.5° | DEFECTO real de coordenadas |
| Cadena libre asistida con ángulos fijos | 120° editorial | 109.5° interno (3 centros medidos) | DEFECTO real de coordenadas |
| Sustituyente auto en anillo rotado 17°, ángulo global fijo 30° | radial/exterior, 120° local | radial error 17°; ángulos locales 137°/103° | DEFECTO real: snap global en sugerencia sin dirección |
| Clic corto en vértice de benceno | un átomo/enlace nuevos | 6 átomos/6 enlaces antes y después | DEFECTO de flujo/hit-test |
| Molécula completa dibujada con arrastres dirigidos alineados | 11 átomos/11 enlaces; enlaces ~40; centros editoriales | 11/11; enlaces 39.99999999999997–40.00000000000004; pares locales observados 119.69–120.32° | PASS para ese gesto; no confirma ni refuta orientación automática |

Imágenes 1×/1.25× del benceno aislado, con sustituyente y molécula completa están en `evidence/before/`. El ciclo cerrado del modelo y la apariencia del renderer se ven regulares; los offsets pi de la representación aromática son deliberados. No se reprodujo defecto de recorte/junction/antialias y no hay motivo basal para tocar `items.py`, QSS o el estilo.

## Decisiones basales

Se aplica la política definida en `design.md`: regla esqueletal planar de 120° separada de la constante tetraédrica 109.5°; click corto automático relativo a geometría local, arrastre explícito conserva snapping y modificadores; cadena zigzag de 120°; anillos existentes se dejan intactos salvo que una regresión numérica futura descubra un problema. No se cambió código antes de capturar este baseline.
