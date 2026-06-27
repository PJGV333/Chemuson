# Cleanup phase 2 architecture report

Rama: `cleanup/architecture-map-and-boundaries`

## Resumen

Segunda fase enfocada en documentación de arquitectura, inventario y diagnóstico. No se hicieron cambios funcionales ni cambios de heurísticas Clean2D/RDKit/ChemName. No se consolidaron más tests históricos porque los restantes contienen casos químicos/regresión no triviales.

## Comandos ejecutados

| Momento | Comando | Resultado |
| --- | --- | --- |
| Baseline | `python -m compileall src tests tools packaging` | Pasa. |
| Baseline | `pytest -q` | 842 passed, 55 skipped, 6 failed en ~55s. Fallos Clean2D conocidos. |
| Baseline | `pytest --collect-only -q` | 903 tests recolectados en ~2.4s. |
| Baseline | `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Pasa. |

## Resultado baseline

El baseline coincide con la fase anterior: pytest falla solo en los 6 tests Clean2D ya conocidos. No apareció ningún fallo nuevo, por lo que la fase continuó.

## Fallos Clean2D que permanecen

| Test | Motivo observado |
| --- | --- |
| `tests/test_clean2d_aromatic_mixed.py::test_naphthalene_fused_system_does_not_collapse` | `run_clean2d_engine` devuelve `complex_preserve_unsafe` en vez de `ok`. |
| `tests/test_clean2d_multilayer_constraints.py::test_tetrandrine_like_hierarchical_blocks_do_not_select_local_graph` | `result.selected` queda en `None`. |
| `tests/test_clean2d_safety.py::test_run_clean_2d_cyclic_geometry_polish_repairs_double_bond_linker` | Longitud del enlace 8-9 queda en ~55.03 frente a esperado ~38.4. |
| `tests/test_clean2d_smart_propose.py::test_smart_clean_repairs_distorted_structure_before_proposing` | Longitud mínima post-clean queda por debajo de `40.0 * 0.70`. |
| `tests/test_clean2d_terminal_rings.py::test_single_anchor_logic_does_not_run_for_fused_rings` | `run_clean2d_engine` devuelve `complex_preserve_unsafe` en vez de `ok`. |
| `tests/test_clean2d_tetrandrine_like_local.py::test_complex_engine_does_not_call_global_redraw_candidates` | El motor reporta `ok` donde el test esperaba preservación/no-op no ok. |

## Inventario de paquetes y archivos

| Archivo / paquete | Tipo | Estado | Recomendación |
| --- | --- | --- | --- |
| `src/chemuson/core/` | producto, API pública | Vivo y de bajo nivel. | Mantener libre de GUI/RDKit. |
| `src/chemuson/chemio/` | producto, API pública e interno | Vivo; contiene conversiones y workers. | Revisar en fase futura dependencias puntuales hacia GUI. |
| `src/chemuson/chemio/_rdkit_worker.py` | interno | Vivo por subprocess. | No borrar por falta de imports entrantes. |
| `src/chemuson/clean2d/` | producto, API pública e interno | Vivo y estratégico; baseline rojo preexistente. | No podar hasta resolver/rebaselinear fallos Clean2D. |
| `src/chemuson/chemcalc/` | producto, API pública | Vivo. | Mantener puro y sin GUI. |
| `src/chemuson/chemname/` | producto, API pública e interno | Vivo; mucha cobertura histórica PR. | Consolidar gradualmente por semántica/acceptance YAML. |
| `src/chemuson/geometry3d/` | producto, API pública e interno | Vivo. | Mantener backends detrás de servicios. |
| `src/chemuson/compchem/` | producto | Vivo. | Mantener sin PyQt; GUI solo en controllers/docks. |
| `src/chemuson/spectroscopy/` | producto, API pública | Vivo MVP. | Mantener predictor desacoplado de docks. |
| `src/chemuson/gui/` | producto, API pública e interno | Vivo; paquete más grande. | Refactor futuro por extracción conservadora desde GUI. |
| `src/chemuson/gui/main_window.py` | producto, GUI orquestadora | Vivo; mezcla ventana, workers, docks y workflows. | Diagnóstico abajo; no partir en esta fase. |
| `src/chemuson/gui/canvas/` | producto, GUI canvas | Vivo; varios archivos grandes. | Refactor futuro por responsabilidades estables. |
| `src/chemuson/gui/controllers/` | producto, controllers GUI | Vivo; buena frontera para orquestación. | Seguir moviendo lógica UI desde `main_window.py` hacia controllers. |
| `src/chemuson/update/` | producto, API pública | Vivo y bien desacoplado. | Mantener sin GUI. |
| `src/chemuson/utils/` | producto/integración | Vivo. | Separar en futuras fases utilidades puras vs integración GUI. |
| `src/chemuson/name2structure/` | producto, API pública | Vivo. | Mantener conectores fuera de GUI. |
| `tools/` | herramienta dev | Vivo; no producto. | No importar desde `src`. |
| `tools/visual/verify_orbitals_rendering.py` | herramienta dev/visual | Vivo opt-in. | Mantener manual o añadir marker CI futuro. |
| `packaging/` | packaging/release | Vivo. | No mezclar con runtime ni pytest normal. |
| `tests/test_clean2d_*.py` | test unitario/integración/regresión | Vivo, crítico. | No consolidar hasta estabilizar baseline. |
| `tests/test_chemname_pr*.py` | test histórico/regresión | Vivo, deuda organizacional. | Clasificación detallada abajo. |
| `tests/data/chemname_acceptance_cases.yml` | acceptance data | Vivo. | Destino natural para casos repetitivos de ChemName. |

Conteos de inventario:

| Grupo | Cantidad |
| --- | --- |
| `src/chemuson/**/*.py` | 172 |
| `src/*.py` | 0 |
| `tests/**/*.py` | 135 |
| `tools/**/*.py` | 6 |
| `packaging/**/*.py` | 8 |

## Clasificación de tests históricos restantes

| Test | Clasificación | Recomendación |
| --- | --- | --- |
| `tests/test_chemname_pr2.py` | Regresión química valiosa de cadena principal. | Renombrar a `test_chemname_core.py` o acceptance YAML si se consolida. |
| `tests/test_chemname_pr3.py` | Regresión de orientación/locantes. | Conservar; candidato a suite `test_chemname_regressions.py`. |
| `tests/test_chemname_pr4.py` | Regresión de sustituyentes halo en alcanos. | Candidato a `test_chemname_functional_groups.py` o YAML. |
| `tests/test_chemname_pr7.py` | Regresión de insaturaciones y alcohol/amina. | Candidato a suite semántica de insaturación/grupos funcionales. |
| `tests/test_chemname_pr8.py` | Regresión anillos/aromaticidad/fusión sin crash. | Conservar; mover con cuidado a suite de anillos. |
| `tests/test_chemname_pr9.py` | Cicloalcanos sustituidos. | Candidato a suite de anillos. |
| `tests/test_chemname_pr10.py` | Ciclohexanos multisustituidos. | Candidato a suite de locantes en anillos. |
| `tests/test_chemname_pr11.py` | Benceno y derivados simples. | Candidato a `test_chemname_special_templates.py` o YAML. |
| `tests/test_chemname_pr12_regression.py` | Regresiones lineales/anillos explícitas. | Conservar como regresión hasta mapear caso por caso. |
| `tests/test_chemname_pr13.py` | Heteroaromáticos simples. | Candidato a suite de heterociclos. |
| `tests/test_chemname_pr14.py` | Diazinas/azoles. | Candidato a suite de heterociclos. |
| `tests/test_chemname_pr15.py` | Naftaleno y sustituidos. | Regresión química valiosa; mover solo con casos intactos. |
| `tests/test_chemname_pr16.py` | Sustituyentes ramificados. | Candidato a grupos funcionales/sustituyentes. |
| `tests/test_chemname_pr17.py` | Fenil/bencil/ciclohexil en cadenas. | Regresión de selección de padre; conservar. |
| `tests/test_chemname_pr18.py` | Selección de padre aromático/cíclico y phenethyl alcohol. | Conservar; candidato a suite de parent selection. |
| `tests/test_chemname_pr19.py` | Naftaleno disustituido. | Regresión valiosa de locantes en fusionados. |
| `tests/test_chemname_pr20.py` | Antraceno/fenantreno. | Regresión valiosa de aromáticos fusionados. |
| `tests/test_chemname_pr21.py` | Quinolina/isoquinolina/indol. | Regresión valiosa de heterofusionados. |
| `tests/test_chemname_pr22.py` | Dienos, triinos, eninos. | Candidato a suite de insaturaciones. |
| `tests/test_chemname_pr23.py` | Aldehído, cetona, ácido, nitrilo. | Candidato a `test_chemname_functional_groups.py`. |
| `tests/test_chemname_pr24.py` | Metoxi, nitro, amino y combinaciones. | Candidato a functional groups/YAML. |
| `tests/test_chemname_pr25.py` | Triazinas/triazoles. | Regresión heterociclos; conservar. |
| `tests/test_chemname_pr26.py` | Benzofurano, benzotiofeno, pireno y locantes. | Regresión valiosa compleja; no consolidar sin revisión química. |
| `tests/test_chemname_pr27.py` | Suite amplia: heterociclos, fused, spiro, bicyclo, estereo. | Dividir en suites semánticas en fase dedicada, no ahora. |
| `tests/test_chemname_pr28.py` | Cargas/isótopos/radicales, sulfoxido, azida, quinonas. | Regresión química valiosa; candidato a `test_chemname_regressions.py` y YAML. |
| `tests/test_aromatic_nitrogen_regression.py` | Regresión GUI/modelo de nitrógeno aromático. | Conservar; no mover hasta separar helpers GUI. |

No se borró ni consolidó ninguno en esta fase porque no hay duplicado trivial evidente sin riesgo de perder cobertura química.

## Diagnóstico GUI y `main_window`

Top 10 archivos GUI más grandes por líneas:

| Archivo | Líneas aprox. | Responsabilidades mezcladas | Extracción futura segura | Extracción peligrosa |
| --- | ---: | --- | --- | --- |
| `src/chemuson/gui/items.py` | 5704 | Items gráficos, geometría de dibujo, interacción visual. | Extraer familias de items por dominio visual con tests snapshot. | Cambiar geometría/wedge/z-order sin baseline visual. |
| `src/chemuson/gui/canvas/canvas_selection.py` | 5154 | Selección, transformación, clipboard, comandos visuales. | Extraer helpers puros de bounds/selección. | Reordenar eventos mouse/undo sin tests UI. |
| `src/chemuson/gui/canvas/canvas_structure.py` | 3875 | Sincronización modelo-escena, import/export visual, análisis, render de estructura. | Separar sync model-scene de análisis químico. | Cambiar serialización o reconstrucción de items. |
| `src/chemuson/gui/orbitals.py` | 3017 | Modelo de orbitales, render iconográfico, paleta. | Separar modelo de specs de funciones de render. | Cambiar coordenadas visuales sin actualizar baseline. |
| `src/chemuson/gui/canvas/canvas_text.py` | 2671 | Texto rico, anotaciones, formato, interacción. | Extraer comandos/helpers de formato. | Cambiar ciclo de edición/foco Qt. |
| `src/chemuson/gui/main_window.py` | 2634 | Ventana, workers, docks, acciones, workflows Name2Structure/CompChem/update. | Mover workers y handlers largos a controllers existentes. | Partir inicialización de docks/actions sin tests de ventana. |
| `src/chemuson/gui/plate_items.py` | 1763 | TLC/gel items, serialización, escala, interacción. | Separar modelos/serialización de QGraphicsItem. | Cambiar bounds/render o escala. |
| `src/chemuson/gui/canvas/canvas_selection_input.py` | 1640 | Eventos de selección mouse/teclado y modo preciso. | Extraer hit-testing puro si ya tiene tests. | Cambiar orden de eventos Qt. |
| `src/chemuson/gui/canvas/canvas_render.py` | 1605 | Export/render de escena, papel, bounds, overlays. | Extraer cálculo de bounds. | Cambiar PNG/SVG/PDF o visibilidad temporal. |
| `src/chemuson/gui/toolbar.py` | 1517 | Toolbars, paletas, iconos, acciones visuales. | Separar builders de paletas por familia. | Cambiar IDs de herramientas o wiring de acciones. |

Recomendación para fase 3 GUI:

| Prioridad | Acción propuesta | Riesgo |
| --- | --- | --- |
| 1 | Extraer de `main_window.py` workers internos (`_DescriptorWorker`, `_NameToStructureWorker`) a módulos controller/service GUI dedicados sin cambiar señales. | Bajo-medio. |
| 2 | Documentar `canvas/` con README corto de mixins y responsabilidades antes de mover código. | Bajo. |
| 3 | Separar helpers puros de cálculo de bounds/selección/render en módulos sin Qt mutable cuando ya haya tests. | Medio. |
| 4 | Dividir `items.py` por familias (`chemical_items`, `annotation_items`, `diagram_items`) solo después de snapshot/visual tests. | Alto si se hace sin baseline visual. |
| 5 | Mantener `Clean2DController` estable hasta resolver los 6 fallos Clean2D. | Alto si se mezcla con fixes. |

## Limpieza suave aplicada

No se aplicaron cambios de código. `ruff` no detectó imports muertos restantes bajo las reglas solicitadas, y cualquier limpieza adicional habría sido cosmética o riesgosa para esta fase.

## Resultado final

| Comando | Resultado final |
| --- | --- |
| `python -m compileall src tests tools packaging` | Pasa. |
| `pytest --collect-only -q` | Pasa; 903 tests recolectados en ~0.8s. |
| `pytest -q` | 842 passed, 55 skipped, 6 failed en ~52s. Misma lista Clean2D que baseline. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Pasa. |

El resultado final conserva exactamente el baseline aceptado: no hay fallos nuevos, no se ocultó ningún test y no se modificó comportamiento funcional.

## Archivos tocados

| Archivo | Motivo |
| --- | --- |
| `docs/architecture.md` | Nuevo mapa de arquitectura y reglas de límites. |
| `docs/cleanup_phase2_architecture_report.md` | Reporte de baseline, inventario, tests históricos y diagnóstico GUI. |

## Riesgos

| Riesgo | Mitigación |
| --- | --- |
| Documentación puede quedar desactualizada si se refactoriza GUI. | Actualizar `docs/architecture.md` en cada fase de extracción. |
| Tests PR de ChemName siguen siendo deuda organizacional. | Migrar por bloques semánticos con tabla de equivalencia caso por caso. |
| Dependencias observadas `chemio -> gui` y `utils -> gui`. | Auditar en fase futura sin cambiar persistencia/autosave de golpe. |

## Próxima fase recomendada

1. Resolver o rebaselinear explícitamente los 6 fallos Clean2D antes de refactors en Clean2D.
2. Crear `docs/gui_canvas_map.md` o sección equivalente para mixins de canvas.
3. Migrar tests ChemName PR a suites semánticas o `tests/data/chemname_acceptance_cases.yml`, conservando cada caso exacto.
4. Auditar dependencias de runtime desde `chemio`/`utils` hacia `gui` y convertirlas en `TYPE_CHECKING` o inversión de control cuando sea seguro.
