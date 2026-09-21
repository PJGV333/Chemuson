# Campaña estratégica de Clean 2D

## 1. Visión

Clean 2D debe convertirse en una característica diferenciadora y confiable de ChemUSON: preservar siempre la química, mejorar de forma medible las estructuras que pueden mejorar con seguridad y degradarse de forma controlada cuando aumenta la complejidad.

Esta campaña es un roadmap y una política de evidencia. No implementa algoritmos. Los contratos actuales de `clean-2d`, calidad, complejidad, snapshots, corpus, baselines y métricas siguen siendo la fuente de comportamiento observable.

La regla maestra es:

- **simple:** una representación simple que ya es correcta no debe empeorar;
- **medium:** las estructuras medianas deben mostrar mejora sostenida y medible a nivel de la familia objetivo y de la campaña;
- **large:** la degradación debe ser controlada, explicable y no destructiva;
- **complex-scale:** la química debe preservarse siempre y debe usarse el mejor candidato seguro entre los candidatos realmente evaluados bajo la estrategia declarada, el espacio de búsqueda y el presupuesto de recursos.

“Mejora sostenida y medible” en `medium` significa mejora reproducible de la familia objetivo y de la campaña; no exige que cada estructura medium mejore en cada cambio. Un caso individual ya bueno puede quedar `unchanged`/`no-op` si no regresa, la familia objetivo muestra mejora agregada o distributiva demostrable y cualquier regresión individual queda visible en el review. Campaign 1 decidirá qué resúmenes son apropiados; esta política no fija todavía un agregado estadístico obligatorio.

“Mejor candidato seguro” no afirma un óptimo matemático ni global: sólo identifica el mejor candidato seguro entre los candidatos efectivamente evaluados bajo los límites declarados.

Si no existe una mejora segura, **preserve-only/no-op** es preferible a un redraw destructivo.

## 2. Principios

1. **La química precede a la estética.** Ningún score compensa una violación de identidad, conectividad, estereoquímica o seguridad.
2. **La evidencia precede a la promoción.** “Se ve mejor” no es suficiente sin baseline, métricas y revisión visual.
3. **El corpus prueba mecanismos generales.** Los casos concretos son regression cases, nunca reglas de ejecución.
4. **La topología precede al pulido.** La optimización local no corrige una decisión global equivocada.
5. **La complejidad exige prudencia.** Un layout seguro y conservador es mejor que una reconstrucción arbitraria.
6. **La observabilidad reutiliza infraestructura existente.** Se extienden `quality_diagnostic`, debug snapshots, baseline reports y baseline diff review; no se crea una segunda telemetría.
7. **Una campaña, un OpenSpec.** No se implementan varias campañas simultáneamente.
8. **Las métricas son vectoriales.** Un scalar puede resumir, pero no ocultar dimensiones ni invalidar hard gates.

### Dirección jerárquica

La dirección arquitectónica futura es:

```text
molecular connectivity
        ↓
topological analysis
        ↓
rigid/semi-rigid block decomposition
        ↓
block graph
        ↓
global block placement
        ↓
connector orientation
        ↓
flexible branch routing
        ↓
internal/local geometry
        ↓
candidate evaluation
        ↓
global ranking
        ↓
local polish
```

La dirección no constituye implementación en esta campaña.

## 3. Las nueve campañas

Cada campaña tendrá su propio OpenSpec con `TARGET FAMILY`, `NON-TARGET FAMILIES`, `EXPECTED IMPROVEMENT`, `HARD INVARIANTS`, `METRICS OBSERVED`, `BASELINE INPUT`, `PROMOTION GATES` y `ROLLBACK CONDITION`.

### Campaign 1 — Benchmark & observability

**Objetivo:** medir antes de modificar.

**Debe producir:** corpus representativo, taxonomía, snapshots opt-in, vector de métricas, baseline report, diff before/after, deterministic replay y clasificación automática de fallos.

**Restricción:** no cambiar el layout productivo, salvo instrumentación estrictamente observacional y explícitamente aprobada.

**Exit criterion:** dos commits pueden compararse y el equipo puede explicar qué cambió, cuánto cambió, en qué familias y con qué decisión.

### Campaign 2 — Topology/decomposition

**Objetivo:** representar una estructura como bloques rígidos/semi-rígidos, conectores y metadata topológica.

**Debe reutilizar:** `MultilayerChemicalGraph`, `BlockGraph`, motifs, `Clean2DComplexityProfile`, `local_graph_cleaner`, `block_unwrap` y los contratos actuales de multilayer constraints.

**Restricción:** no crear un modelo paralelo si `core` ya contiene información equivalente. Primero se validan tests topológicos; la geometría permanece igual durante la fase inicial.

**Exit criterion:** una estructura medium/large puede descomponerse determinísticamente en bloques, conectores y metadata.

### Campaign 3 — Medium molecule assembly

**Objetivo:** mejorar estructuras de aproximadamente 20–60 átomos del grafo con múltiples anillos o cadenas.

El global placement precede al local polish. La evidencia debe mostrar mejora sostenida y reproducible de la familia medium objetivo y de la campaña, sin regresión de los casos simples. Un caso medium individual ya correcto puede permanecer unchanged/no-op si su resultado no regresa y la evidencia agregada o distributiva de la familia queda documentada; cualquier regresión individual debe ser visible en el review.

**Exit criterion:** mejora reproducible en la familia medium y Gate B satisfecha.

### Campaign 4 — Rigid/fused/multiring layout

**Objetivo:** orientar bloques rígidos, sistemas fused, spiro, bridged, multiring y sustitución congestionada preservando la geometría interna cuando sea razonable.

La campaña debe distinguir claramente naphthalene-like, fused, spiro y bridge como familias topológicas, no como nombres de moléculas.

### Campaign 5 — Flexible connectors and branch routing

**Objetivo:** separar ramas, evitar cruces, orientar conectores, distribuir ocupación angular e impedir que fragmentos grandes compitan por el mismo espacio.

Las decisiones de routing deberán ser expresables por propiedades de conectores y del block graph; no por IDs de fixtures.

### Campaign 6 — Macrocycles and large structures

**Objetivo:** layouts jerárquicos para macrocycles, estructuras large multiblock, natural-product-like, peptide-like y glycoside-like.

No se permitirá global redraw destructivo si no existe una alternativa segura y aprobada. La política compleja existente y `preserve-only` mantienen prioridad.

### Campaign 7 — Global candidate search/ranking

**Objetivo:** generar varias disposiciones topológicamente razonables.

El orden es obligatorio: hard gates primero y soft metrics después. El proceso debe permanecer determinista, limitar el crecimiento combinatorio y conservar candidate source, candidate count y motivo de selección.

La disponibilidad o éxito de RDKit, CoordGen u otro backend no define la calidad; solo aporta candidatos.

### Campaign 8 — Local polish

**Objetivo:** refinar bond lengths, angles, minor overlaps, label separation y geometría local sin reorganizar la decisión topológica global.

`local polish SHALL NOT` reorganizar bloques globales salvo contrato explícito. Todo cambio debe preservar stereo layout, límites de desplazamiento y contratos de selección.

### Campaign 9 — Production acceptance

**Objetivo:** validar corpus amplio, determinismo, performance, estabilidad de regresión y comparación visual manual.

Debe definirse una versión identificable del baseline de producción.

**Exit criterion:** Clean 2D puede considerarse una capacidad diferenciadora y confiable de ChemUSON con rollback documentado.

## 4. Contratos de seguridad

### Hard constraints

Un candidato que viola cualquiera de estas propiedades no puede ganar por tener mejor estética:

- atom IDs preservados;
- bond IDs preservados;
- atom count y bond count preservados;
- element identity y formal charge preservadas;
- bond endpoints y bond order preservados;
- aromaticity preservada;
- stereo metadata preservada;
- selection metadata preservada cuando aplique;
- coordenadas completas y finitas;
- integridad del `MolGraph` y de sus connected components.

Se reutilizan los estados existentes `applied`, `rejected`, `preserve-only`, `no-op` y `failed-controlled`, y las razones estables existentes: `invalid-coordinates`, `invariant-violation`, `stereo-risk`, `boundary-bond-risk`, `new-crossing-risk`, `collision-risk`, `collapsed-ring-risk`, `excessive-displacement`, `worse-quality` y `backend-failure`.

### Hard gates

Los hard gates se ejecutan antes del ranking:

1. invariantes del grafo;
2. coordenadas válidas y finitas;
3. preservación de identidad química;
4. preservación estereoquímica;
5. integridad de bonds de frontera;
6. ausencia de cruces estructurales nuevos;
7. ausencia de colisiones inseguras;
8. rings no degenerados;
9. desplazamiento dentro del contrato de la estrategia;
10. resultado controlado ante fallo de backend.

### Soft metrics

Solo los candidatos que sobreviven se comparan por métricas. El vector puede incluir:

- `bond_length_error`;
- `bond_length_variance`;
- `bond_angle_penalty`;
- `atom_collision_count`;
- `label_collision_count`;
- `bond_crossing_count`;
- `ring_distortion`;
- `ring_degeneracy`;
- `rigid_block_distortion`;
- `branch_separation`;
- `connector_congestion`;
- `compactness`;
- `whitespace_balance`;
- `global_extent`;
- `symmetry_preservation`;
- `candidate_count`;
- `candidate_source`;
- `runtime_ms`.

Las definiciones, unidades, polaridad, optionality y tolerancias ya existentes son autoritativas. No se fijan thresholds numéricos nuevos sin baseline suficiente.

## 5. Taxonomía de moléculas

La clasificación es ortogonal. Cada caso tiene una `size class` y tags de familia/topología.

### Size class

- `simple`;
- `medium`;
- `large`;
- `complex-scale`.

`medium` tiene como objetivo inicial aproximadamente 20–60 átomos del grafo. Atom count nunca es el único criterio.

### Topology/family tags

La taxonomía extensible contempla:

`acyclic`, `branched`, `monocycle`, `aromatic`, `multisubstituted-aromatic`, `fused`, `spiro`, `bridged`, `macrocycle`, `multiblock`, `rigid-flexible`, `peptide-like`, `glycoside-like`, `heteroatom-rich`, `charged`, `stereo-sensitive`, `coordination`, `selection-boundary`, `congested`, `known-delicate`, `known-failure`.

Las etiquetas actuales como `baseline`, `known_delicate`, `complex_policy_guard` y `stereo_sensitive` se conservan para compatibilidad.

### Metadata medible

Cuando esté disponible se registra:

`atom_count`, `heavy_atom_count`, `bond_count`, `ring_count`, `connected_components`, `rigid_block_count`, `rotatable_connector_count` y `macrocycle_count`.

La ausencia inicial de un campo no cambia el significado de un case ID.

## 6. Identidad del corpus

Cada regression/benchmark case tiene un ID estable, único y apto para pytest IDs y snapshots. La identidad del caso se separa de `current expected quality`.

No se renombra un `known_failure` para ocultar que cambió. Cambios en fixture, tags, modo, target o expected states aparecen explícitamente en el diff. Los casos `known_delicate`, `selection_boundary`, `stereo_sensitive` y `complex_policy_guard` permanecen visibles y observacionales.

El corpus actual se usa como evidencia de mecanismos: simple, aromatic, fused, heteroaromatic, charged, stereo-sensitive, selection-boundary, multi-block y macrocycle. Las cifras de cobertura se obtienen del registro actual, no se codifican en producción.

## 7. Métricas

Las métricas responden a dos preguntas separadas:

- **Seguridad:** ¿pasa los hard gates?
- **Calidad:** entre los candidatos seguros, ¿qué cambió y por cuánto?

Las métricas de geometría son diagnostic-only hasta que un OpenSpec posterior promueva explícitamente una métrica. `visual_score` conserva su semántica existente de lower-is-better cuando se use el score interno; el score normalizado de `quality_reporting` es un valor de reporting 0–1 y no altera ranking.

Valores ausentes solo se permiten cuando la definición declara que el caso es no aplicable o unavailable. No se aceptan NaN ni infinity en reports JSON.

## 8. Benchmark y baseline policy

Antes de cualquier cambio algorítmico:

1. capturar baseline;
2. ejecutar el corpus idéntico;
3. persistir el report;
4. modificar el algoritmo;
5. ejecutar el corpus idéntico;
6. producir diff before/after;
7. clasificar regresiones;
8. revisar visualmente casos relevantes;
9. aceptar o rechazar el cambio.

El baseline no se actualiza para hacer verde un test. Una actualización necesita justificación explícita, identidad de corpus, motivo, campos afectados y aprobación de la revisión.

La comparación debe responder: **What improved? What regressed? By how much? On which families? Why was the candidate selected? Which hard gates were checked?**

## 9. Promotion gates

Una promoción importante requiere:

- **Gate A — chemical safety:** cero violaciones nuevas de invariantes;
- **Gate B — simple non-regression:** los casos simples correctos no empeoran;
- **Gate C — target-family improvement:** la familia objetivo mejora reproduciblemente;
- **Gate D — cross-family regression review:** toda regresión no objetivo queda identificada y justificada;
- **Gate E — determinism:** no aparece aleatoriedad no controlada;
- **Gate F — performance:** no existe degradación injustificada ni crecimiento combinatorio descontrolado;
- **Gate G — manual visual review:** casos before/after representativos revisados.

No se fijan porcentajes arbitrarios en Campaign Policy. Los porcentajes, si algún día resultan necesarios, deben venir de evidencia y otro OpenSpec.

## 10. Protocolo experimental

Una estrategia experimental declara su routing, candidate source, parámetros, seed, familias objetivo y no objetivo, invariantes, métricas y condición de rollback. Puede ejecutarse en benchmark o detrás de routing explícito; no se vuelve default automáticamente.

Los backends externos son candidate sources. `backend success != layout accepted`. Toda salida se evalúa con contratos propios de ChemUSON.

La revisión visual humana complementa las métricas; no las reemplaza. Debe usar casos representativos, misma escala cuando proceda, anotaciones de cambios relevantes y vínculo con el baseline report.

## 11. Política de determinismo

Para las mismas entradas, modo, parámetros y seed, el resultado y el orden de candidatos deben ser reproducibles dentro de tolerancias documentadas. Si un backend no es determinista, se registra backend, seed si existe, source, candidate ordering observado y decisión final auditada.

Debug snapshots siguen siendo opt-in mediante la infraestructura existente. Desactivarlos no puede cambiar routing, ranking, selección ni geometría.

## 12. Política de performance

Campaign 1 mide wall time, candidate generation time y candidate count por size class y topología. Se separan tiempos de backend, evaluación y ranking cuando sea posible. El objetivo inicial es observar distribución y crecimiento, no imponer un SLA inventado.

Una estrategia que genere crecimiento combinatorio sin control no se promociona. Cualquier límite futuro debe justificarse por corpus y evidencia, no por una excepción a un caso concreto.

## 13. Protocolo de revisión visual

La revisión visual acompaña al diff cuantitativo en cambios importantes:

1. seleccionar casos objetivo y no objetivo representativos;
2. renderizar before/after con parámetros documentados;
3. comprobar identidad, selección, stereo y boundary bonds;
4. observar crossings, collisions, ring shape, branches, whitespace y labels;
5. anotar mejoras y regresiones con case ID estable;
6. relacionar cada observación con una métrica o declararla como juicio humano separado;
7. conservar la decisión en el reporte de la campaña.

“Se ve mejor” sin baseline, métrica o registro visual no es criterio de promoción.

## 14. Rollback policy

Se revierte o se mantiene fuera de production routing una estrategia que:

- viole química o hard constraints;
- introduzca regresiones en casos simples correctos;
- empeore ampliamente familias no objetivo;
- pierda determinismo;
- cause degradación severa e injustificada de performance.

No se responde primero agregando excepciones por molécula. Se conserva la evidencia como test, report o experimento cuando sea útil para una estrategia posterior.

## 15. Definición de production ready

Una campaña está lista para producción solo cuando:

- su OpenSpec individual está aprobado y validado;
- el corpus y la taxonomía son reproducibles;
- la identidad de casos y baseline version están identificadas;
- todos los hard constraints y Gate A pasan;
- los casos simples pasan Gate B;
- la familia objetivo pasa Gate C;
- las regresiones no objetivo están revisadas;
- determinismo y performance están medidos;
- la revisión visual manual está registrada;
- existe una condición de rollback practicable;
- los reportes, diagnósticos y fuentes de candidatos explican la decisión.

Production ready no significa que toda molécula tenga una representación ideal. Significa que el sistema preserva la química, mejora donde tiene evidencia de hacerlo con seguridad y degrada de manera controlada y auditable.
