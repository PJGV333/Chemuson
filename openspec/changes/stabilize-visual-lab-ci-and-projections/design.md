# Diseño — CI por duración y proyecciones con semántica explícita

## Evidencia de baseline

Rama `feat/molecular-depiction-visual-lab`, SHA inicial `87e298c46535e1dffb0477f36412b3dedc87bd98`, árbol limpio y sincronizado con origin. El baseline obligatorio está en `baseline-command-output.txt`; el plan y artefactos de Actions en `baseline/ci-artifacts/`. El análisis completo de los dos timeouts y los tiempos por archivo está en `baseline.md`. OpenSpec estricto debe pasar antes de tocar producción o herramientas.

## Sharding

Usar Largest Processing Time (LPT) sobre unidades indivisibles por archivo. El peso es `node_count × seconds_per_test` de `tests/ci/pytest_file_durations.json`; archivos sin historia reciben un fallback explícito y versionado. El orden de desempate usa peso, cantidad, primer ordinal de colección, ruta y shard index, todos deterministas. Preservar el orden original de node IDs dentro de cada archivo y garantizar asignación exacta una vez.

El plan registra la huella del perfil y la duración estimada por shard. El runner habilita un plugin local que escribe una línea por test terminado y agrega segundos por archivo al JSON de resultado, incluso cuando otro test posterior agota el timeout. CI sube también el JSONL. La verificación sigue siendo fail-closed para timeout, código de salida, asignación y cobertura; el límite externo permanece en 300 s y la limpieza POSIX sigue terminando el grupo completo.

## Viewport aislado del laboratorio

`ProductionCanvasAdapter` usa una subclase interna del canvas destinada al laboratorio. La subclase sobreescribe sólo los puntos de centrado inicial diferido y resize para centrar `graph_scene_rect`, no `paper_width/paper_height`; la ruta normal de `ChemusonCanvas` queda intacta. Se aplica el mismo transform, centro de escena y sincronización de scrollbars a Current y Working después de construir casos, zoom/rotación y resize. No se usa monkey-patching ni se añade estado al grafo.

Las pruebas mostrarán una ventana Qt real, esperarán los timers/eventos y comprobarán pixeles no vacíos, centros de ítems dentro del viewport, conteo de átomos/enlaces visibles, geometría/escala coincidentes y repetición tras resize/cambio de caso/zoom/rotación.

## Geometría Haworth

Definir coordenadas canónicas no regulares independientes para anillo de seis y cinco miembros; numeración estable desde el oxígeno y aristas frontales por índices canónicos, no por umbral sobre coordenadas transformadas. Aplicar aspecto, perspectiva, skew y rotación sobre el modelo canónico con rangos finitos/documentados. Una tabla específica `(ring_size, ring_position, up/down)` determina vectores de los grupos, despeja el exterior y corrige especialmente C3 OH arriba. Los defaults beta-D-glucopiranosa conservan C1 OH arriba, C2 abajo, C3 arriba, C4 abajo y C5 CH2OH arriba. El renderer de producción validará colisiones de label/path a varias escalas; AB3 queda sin cambios.

## Contrato de estereoquímica

Añadir tipos inmutables de descriptores de proyección y un `ProjectionInstance` que contenga la geometría, un `MolGraph` fresco de sólo representación y los descriptores tipados. Los descriptores guardan índices de átomos/posiciones y no se convierten en `BondStereo`/`stereo_ez` inventados. Haworth registra α/β, D/L y orientación de cada sustituyente; Fischer registra la convención horizontal-hacia-observador y vertical-alejándose, lados y centro indexado, sin inferir R/S cuando no se declara; silla registra conformación y para cada grupo su relación up/down y axial/equatorial.

`MolGraph`/`.cmsn` no serializan el contenedor de proyección. Por tanto, la instancia tipada es el contrato seguro para conservar semántica en memoria; la adaptación legacy `build_*_template() -> MolGraph` y `ProjectionGeometry.to_molgraph()` son sólo dibujos y no prometen reconstruir estereoquímica. No se toca el formato `.cmsn`; la pérdida de descriptor sólo puede ocurrir al escoger explícitamente el adaptador gráfico, y queda documentada y probada.
