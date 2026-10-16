# Diseño — Molecular Depiction Visual Lab

## Baseline y alcance

La rama de campaña parte de `origin/feat/molecular-depiction-excellence` en `ee5294800c5779f7451e625d8089025e8e6632ab`. El árbol estaba limpio; el remoto experimental fue actualizado mediante fetch. El inventario exacto de comandos está en `baseline-command-output.txt`; resumen en `baseline.md`. El renderer a calibrar es el estado de esa base, no los cambios experimentales del laboratorio.

## Renderer compartido

La interfaz aloja dos `ChemusonCanvas` reales. Cada caso proporciona un `MolGraph` inmutable de referencia y ambos adaptadores crean copias profundas independientes del mismo grafo; los canvas ejecutan `_rebuild_items_from_model`, el cálculo de contexto de anillos/joins y el trimming de etiquetas que usa producción. La vista original conserva `CHEMDOODLE_LIKE`; la vista ajustada aplica un `DrawingStyle` derivado del perfil Working a sus `AtomItem`/`BondItem` existentes y recalcula los trims con los métodos de canvas existentes. No hay un segundo algoritmo de dibujo.

Los datos experimentales no se aplican a `ChemusonWindow`, QSettings del usuario, documentos, `MolGraph` o UndoStack de producción. Se limita la extensión de renderer a campos opcionales en `DrawingStyle`, cuyo valor `None` deja intactos los cálculos actuales: clearance de etiqueta normalizado, anchura relativa de punta de cuña, desplazamiento terminal hash para encuentro sigma, espaciado explícito de dobles bold y espaciado de triples. Se mantiene `BondItem` API pública actual y los defaults visuales.

## Modelo y controles

`DepictionStyleParameters` guarda proporciones relativas a la longitud de enlace de referencia de `CHEMDOODLE_LIKE` (40 scene units), el recuento hash como entero y `cap_style`/`join_style` como opciones discretas. Parámetros directos soportados: grosor de trazo base; offset de doble enlace e inner-trim; anchura de cuña; cantidad, ancho mínimo/máximo y grosor de hash; opciones discretas cap/join. Controles experimentales opt-in: clearance etiqueta/enlace, tip de cuña, desplazamiento terminal hash, offset de doble bold y offset triple. Los overrides opt-in se guardan como `null` en Current (significa política del renderer actual); al activarlos la GUI muestra y modifica el valor numérico real.

Los límites son acotados para evitar geometría no finita, offsets negativos o polígono degenerado. El validador rechaza campos desconocidos y tipos ambiguos. JSON usa `schema=chemuson.depiction-style`, `version=1`; migraciones futuras se registran explícitamente por versión y perfiles futuros no reconocidos fallan cerrados. La ruta de perfiles vive en la ubicación de configuración propia del laboratorio, nunca en preferencias normales de ChemUSON.

No se exponen como controles editables los factores internos de miter/fork de cuña, la anchura estrecha calculada del tip con vecinos, orientación de la línea pi, tamaño de fuente/side bearings, trim geométrico de rings, ni caps per-style: actualmente dependen de topología, etiquetas o constantes internas y no se pueden aislar de forma segura. Se listan como no disponibles en documentación/UI, no como sliders falsos.

## Casos y orientación

El corpus tendrá 12+ IDs estables (A1–A4, B1–B4, C1–C4). Las moléculas con valencia y estereoquímica coherentes se etiquetan como moléculas/fixtures de depiction; los escenarios geométricos abstractos se marcan explícitamente como fixtures geométricos. La escena Haworth se construye con la plantilla determinista existente. Cada render usa dos grafos clonados, y pruebas comparan sus snapshots antes/después. La orientación de estrés usa la transformación de vista Qt compartida (no toca coordenadas ni `BondStereo`); no se transforma químicamente el grafo.

## Interfaz, diagnóstico y exportación

La ventana muestra dos paneles de igual tamaño, selección de caso, tema, zoom 50–300 %, ángulo de vista, controles agrupados, y botones no modales para Current/Working, reset, perfiles, diagnóstico y exportación. Scroll/zoom/rotación se sincronizan. Los overlays de diagnóstico se crean como ítems de escena descartables: centros atómicos, extremos visibles, paths/contornos de cuñas, cajas de label y longitudes/gaps. Los exportadores construyen escenas sin overlays, salvo petición explícita.

La exportación estable incluye PNG Current, PNG Working, comparación, hoja de contacto del corpus, perfil JSON y Markdown con IDs, parámetros, versión y SHA. Render utiliza tamaños, fuente elegida, fondo, encuadre y antialiasing registrados. SVG sólo se añadirá si puede renderizarse mediante la escena Qt compartida sin duplicar el motor; de otro modo se documenta como no disponible.

## Integridad arquitectónica y promoción futura

El laboratorio vive bajo `tools/visual`, no es un paquete runtime `src/chemuson` ni importa nuevos subsistemas de dominio. Depende de GUI existente desde herramienta de desarrollo; no se introduce una dependencia inversa ni circular. El siguiente paso de promoción, tras aprobación visual, será una propuesta OpenSpec separada que traslade sólo los valores aprobados a una configuración compartida de `DrawingStyle`, conserve los defaults/compatibilidad solicitados y cubra exportación/renderer. El laboratorio no modifica valores de producción por sí mismo.

## Riesgos y límites

- `ChemusonCanvas` es una vista de edición grande y su constructor crea overlays; la herramienta debe desactivar interacción/papel/grid y cerrar explícitamente canvas/ventana en pruebas.
- Los defaults de color tipográfico están definidos por el renderer; tema oscuro requiere adaptar colores de presentación de ambos paneles por igual sin incluir esos colores en la calibración geométrica.
- Las capturas pueden variar con Qt, fuentes o DPI. Tests priorizan geometría, invariantes y hashes dentro del mismo entorno; los controles visuales usan tolerancias y reporte de entorno, sin snapshots de píxel frágiles.
- El conjunto de pruebas CI histórico de la rama experimental documentó fallo dependiente del orden en `test_compchem_controller_generates_async_with_fake_backend` en shard 4 y shards 5–8 no ejecutados. Este estado no se atribuye al laboratorio; validar mediante los shards actuales y reportar los resultados exactos.
