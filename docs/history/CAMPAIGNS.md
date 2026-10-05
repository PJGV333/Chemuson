# Memoria de campañas de ChemUSON

**Propósito:** conservar decisiones, resultados y aprendizaje reutilizable de campañas cerradas sin usar ramas, informes de ejecución o capturas intermedias como museo. Los contratos OpenSpec archivados conservan los requisitos detallados; este archivo es el índice narrativo. No autoriza trabajo nuevo ni sustituye un OpenSpec activo.

## Punto de referencia

- Fuente examinada: `main` y `origin/main` en `1db4f63b52af79247745b3a8a220fb728348218c` (QA de modernización UI, 2026-10-03), antes de esta campaña de higiene.
- La rama de higiene se crea desde ese SHA; no se modifica ni se integra a `main` en este trabajo.
- Versión del producto: `0.3.0-dev`, sin bump ni publicación.
- Las listas completas de ramas/SHAs antes y después de la poda, las métricas y la verificación de higiene están en [`REPOSITORY_CLEANUP_2026-10-03.md`](REPOSITORY_CLEANUP_2026-10-03.md).

## Campañas de arquitectura y mantenimiento

### Límites de módulos, resiliencia y editor — julio a septiembre de 2026

Las campañas OpenSpec archivadas construyeron el catálogo de módulos y movieron responsabilidades con cambios incrementales: desacople de autosave/resiliencia, vista molecular compartida, límites de persistencia/GUI, raíz de composición, shell de aplicación, workers, geometría Clean2D de ventana y selección del canvas. Las fases posteriores documentaron M21 (`platform.settings`), M22 (`resilience`), M19 (composition root) y M20 (selección de `editor2d`). Se evitaron renombrados masivos y se mantuvieron shims solo cuando había consumidores de compatibilidad.

**Cierre observado en la campaña de selección:** suite de arquitectura hasta 268 tests, suite completa 1496 passed/55 skipped, compileall y Ruff scoped correctos. Un `openspec validate` anterior tenía un fallo de requisito en `application-composition-root`; fue corregido en el cierre posterior. El baseline siguió evolucionando: esos números históricos no describen la rama actual.

**Lecciones:** auditar consumidores reales antes de mover/eliminar un módulo; distinguir fachada compatible de lógica activa; conservar el orden de despacho Qt y el ownership de selección; actualizar `architecture/modules.yml` junto con dependencias/paquetes; nunca usar un nuevo baseline para ocultar una regresión.

### Limpieza conservadora anterior

- `cleanup/remove-dead-code-and-tests` — punta integrada `933a27cb068b22f251fdad354e4c63cde8e875ee`: retiró stubs/scripts de depuración sin consumidores, imports muertos y wrappers de tests redundantes cuando existía cobertura semántica equivalente; conservó los casos químicos no triviales. La utilidad visual orbital se mantuvo como herramienta opt-in (`tools/visual/verify_orbitals_rendering.py`). En aquel baseline se documentaron seis fallos Clean2D preexistentes; no representan el estado actual.
- `cleanup/architecture-map-and-boundaries` — `71798871e4298b24130a6502ed458dd512744c6c`: inventario y mapa, sin cambio funcional; recomendó no consolidar tests ChemName/ Clean2D cuando sus casos codifican regresiones químicas distintas.

El detalle voluminoso de esos reportes quedó consolidado aquí; los contratos y OpenSpecs archivados siguen siendo la especificación, no se eliminaron tests químicos actuales.

## Modernización UI — Fases 0–8

**Cierre integrado:** `1db4f63b52af79247745b3a8a220fb728348218c`; rama de QA histórica `release/ui-modernization-qa` terminaba en el mismo SHA. Siete cambios OpenSpec están archivados bajo `openspec/changes/archive/2026-10-03-*`; la Fase 8 documenta QA y permanece como cierre activo. El protocolo estricto validó 42 cambios, sin fallos.

| Fase | Resultado / decisión conservada |
|---|---|
| 0–1 — fundación visual | Se aprobó un prototipo PyQt6 como contrato visual; se centralizaron tokens, temas claro/oscuro, preferencias y acceso a iconos. `theme.py` se conserva porque la especificación global lo identifica como referencia exacta de diseño. |
| 2 — iconos | Migración SVG/HiDPI usando infraestructura existente y recursos incluidos; sin dependencia externa. |
| 3–4 — rail y convergencia | `QWidget` rail (58 px), botones de 42 px, 15 acciones en seis grupos; `QMenuBar` y toolbars clásicas ocultas, no retiradas; flyouts reutilizan `QAction` vivos. Menú Alt/hamburguesa, flyout en segundo clic/clic derecho. La ventana encaja a 980×600. Comparación numérica 13/13 y matriz visual 11/11. |
| 5 — panel y estado | Siete docks existentes se presentan en un SidePanel; se conservan widgets/señales/handlers. Estado de 34 px; páginas visibles y overflow documentados. |
| 6 — app bar y documentos | Barra de aplicación y pestañas de documentos; se conserva edición/undo/redo y navegación existente. |
| 7 — CommandPalette y polish | Paleta `Ctrl+P`; `Ctrl+K` sigue reservado a Clean2D. Se verificaron onboarding, plantillas, preferencias, flyouts y HiDPI; selección de plantillas con clic simple. |
| 8 — QA y cierre | Manual, notas de release, cinco capturas finales, Flatpak instalado y round-trip `.cmsn` aprobados. PyInstaller arranca. El script existente llama `.AppImage` a un ejecutable portable; sin `appimagetool` no se afirmó un contenedor AppImage Type 2. |

**Evidencia QA registrada en el cierre:** compileall pasó; 1816 tests recolectados; suite completa 1760 passed, 55 skipped y un fallo conocido (`tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`), idéntico a baseline; arquitectura 269 passed; UI dirigida 304 passed; Flatpak smoke PASS. Ruff scoped mantiene un F401 no relacionado en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`. El aviso de `propagateSizeHints()` bajo el plugin Qt offscreen no fue fatal.

**Decisiones/no objetivos:** sin cambios a lógica química, Clean2D, ChemName, geometría de plantillas, serialización `.cmsn`, versión o dependencias; sin iniciar Fase 9. Las capturas de `docs/ui-modernization/after/` son la evidencia canónica de estado final; las capturas reales aprobadas KDE/Wayland se conservan aparte. Iteraciones intermedias, generadores de una sola fase y el ejecutable del spike se retiraron tras resumir decisiones y resultados.

### AI Molecular Structure Bridge — Phase 1 Foundation y Phase 2 UI mínima

Phase 1 se cerró y archivó como `2026-10-04-define-ai-molecular-structure-bridge`. M23 ofrece servicio provider-neutral, salida JSON SMILES estricta, validación ChemIO aislada y adaptador HTTP OpenAI-compatible explícito; no modificó Clean2D. Phase 2 se implementa en `add-ai-molecular-assistant-ui`: QAction en Structure/Command Palette, formulario modeless con configuración explícita/transitoria, controlador M10 con worker en hilo aparte, revisión del resultado y confirmación antes de inserción mediante el macro undoable existente. M10 consume M23; M23 mantiene únicamente M00/M01 y Clean2D sigue independiente. Las pruebas offline verifican fallos/abandono sin mutación, ruta de inserción y undo/redo normal. Por el límite solicitado de 5–10 minutos, no se repitió la suite completa: el baseline tardó 19:26, con el fallo histórico `test_compchem_controller_generates_async_with_fake_backend`; Ruff mantiene el F401 histórico de la prueba Clean2D. Phases 3–5 siguen sin iniciar en este punto.

### Clean2D — campañas separadas y no integradas

**No se modificó código Clean2D en la campaña UI ni en esta higiene.** En `main` se conservan las campañas OpenSpec de julio sobre corpus, snapshots, métricas, baseline/diff review, determinismo, preservación compleja y layouts aromáticos/mistos. Sus propuestas y baselines permanecen en los archivos OpenSpec archivados.

Hay trabajo posterior valioso que no forma parte de `main` y no debe tratarse como integrado:

| Rama remota retenida | SHA observado | Trabajo único observado | Decisión |
|---|---|---|---|
| `origin/clean2d/campaign-implementation` | `6b153e3699bf2b04e60189a4606040f766c8eb5c` | 18 commits únicos; evidencia/código/tests para campañas 1–4 (observabilidad, descomposición topológica, ensamblaje medium y multiring). | Retener sin integrar ni borrar; revisar cada campaña/OpenSpec antes de promoción. |
| `origin/clean2d/campaign-policy` | `31c565c8dfe15f2b4c2287eb4ece56e2b67d9732` | Dos commits únicos de política. | Retener; la política es contexto, no autorización para aplicar algoritmos. |
| `origin/ui/modernization` | `677b8294138357774e3937ceecbb97dd621b8f26` | Historia de la UI más trabajo Clean2D no integrado; contiene `docs/clean2d/CAMPAIGN.md` y cambios de contrato/corpus. | Retener; UI integrada no convierte los commits Clean2D exclusivos en integrados. |

La dirección estratégica del documento `docs/clean2d/CAMPAIGN.md` de esa historia queda resumida así: medir antes de modificar; preservar identidad, conectividad y estereoquímica como hard gates; descomponer topología antes del pulido local; comparar candidatos seguros con métricas vectoriales; mantener `preserve-only/no-op` ante falta de mejora segura; establecer cada campaña con familia objetivo/no objetivo, baseline, gates, determinismo, coste y rollback. Las nueve campañas planificadas abarcan observabilidad, descomposición, medium, multiring, conectores, macrociclos, búsqueda/ranking, pulido local y aceptación de producción. Esta síntesis no incorpora aquella rama ni altera los contratos vigentes.

**Experimentos de parche no adoptados:** los archivos raíz `clean2d_failed_local_graph_attempt.patch` y `clean2d_failed_tetrandrine_selection_integrity.patch` son parches, no módulos importados. Documentan intentos tempranos de enrutar un limpiador local de grafo y endurecer la integridad de selección/estereo en casos complejos; no se aplican automáticamente. Se retiran como artefactos sueltos, preservando su contexto en este resumen y los objetos Git históricos (commit `aeaef0585bebf4f9bf64b4afdcc1e4c90a3f6753`, 2026-06-18). La implementación actual `src/chemuson/clean2d/local_graph_cleaner.py` no se toca.

## Otros trabajos únicos que se mantienen fuera de `main`

No se borran ramas no-ancestro por antigüedad o nombre. En el reporte de cierre consta inventario completo con SHAs; requieren revisión del propietario antes de cualquier decisión:

- `origin/docs/chemuson-interactive-article` (`cdcae3c2e37df6967c61c7e7b315db0e225cfb67`): 36 commits únicos; artículo interactivo y assets, pendiente decidir publicación/integración.
- `origin/refactor/packaging-pyproject` (`d44115b349fdd4be3bd25599c89cf605cb33613b`): dos commits únicos de runtime hook PyInstaller; verificar frente al empaquetado actual.
- `origin/gh-pages` (`c856fc1445a11576a5a2718cd40380bc54056853`): publica el remoto Flatpak. Su árbol actual contiene 6895 objetos `flatpak/beta/repo/objects/*` (159,512,305 bytes lógicos); no eliminarlos ni reescribirlos en una tarea de higiene.
- `origin/codex/implementar-sistema-de-apariencia-en-chemuson` (`2a40df84137ca07bc4fd2583fd3ec31ff9b0eb47`): commit único de tema/apariencia, requiere comparación antes de descartar.
- `origin/tema-chemuson` (`167ade37081140747bb014cc613646ffdcf3063d`): cuatro commits únicos de trabajo visual; conservar hasta cotejar decisiones.

## Fallos y decisiones conservadas

- Un nombre de rama o una antigüedad no demuestra que sea prescindible: se comparan SHAs, ancestros y árboles. Las ramas únicas se conservan hasta que su autoría/propósito se revise; las ramas redundantes solo se eliminan después de comprobar ancestro de `origin/main` y registrar los SHAs.
- El prototipo visual y sus valores aprobados guiaron producción; scripts/código de demostración sin consumidor se retiraron después de conservar los tokens normativos, capturas útiles y lecciones Qt. `Ctrl+P`, menús ocultos, orden de despacho, `QSettings`, viewport de `QScrollArea` y DPR se validan mediante contratos/tests, no solo imágenes.
- El archive visual orbital era salida generada de comparación, no fixture: usaba una referencia absoluta de otra máquina y múltiples subreportes. El resumen principal daba `average_ref_iou=0.6461`; el baseline de una comparación quirúrgica daba `average IoU=0.7063`, con métricas distintas por familia. No debe interpretarse como baseline de producto vigente ni como autorización para cambiar geometría orbital. Para no recrear artefactos en `tests/data`, `tools/orbital_fit_report.py` y `tools/render_orbital_family_preview.py` escriben por defecto en el directorio temporal del sistema (`chemuson/orbital-fit-report`); el primero acepta `--output-dir`, ambos aceptan `CHEMUSON_ORBITAL_REPORT_DIR` y el preview también conserva sus flags de salida explícita. El preview pasó el smoke externo; el reporte alcanza un fallo preexistente de `PiBondingParams.ring` en `_family_metric_strings`, tras generar salidas parciales solo en `/tmp`. Se conserva como `NEEDS_REVIEW` por su posible uso manual, sin reparar una cuestión orbital fuera de alcance.
- `src/sys` era PostScript de ImageMagick de 11,708,416 bytes sin consumidores/imports/entrypoints/packaging; `src/repro_v2.png` era un PNG de reproducción sin referencias. Se retiran del árbol actual, sin tocar ningún algoritmo químico. Git history no se reescribe.

## Pendientes explícitamente no iniciados

- Resolver el único fallo de baseline de CompChem y el F401 de la prueba Clean2D requieren tareas OpenSpec separadas; no se arreglan aquí.
- Comparar o integrar ramas únicas y resolver el almacenamiento histórico del remoto Flatpak requiere decisión del propietario.
- Cualquier nueva fase UI, modificación de Clean2D/ChemName, bump de versión, publicación o reescritura de historia necesita alcance y aprobación propios.
