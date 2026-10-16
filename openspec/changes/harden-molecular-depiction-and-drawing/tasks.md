# Tareas — robustez de depiction y dibujo

## 1. Rama y baseline
- [x] 1.1 Confirmar rama experimental/upstream, SHA inicial y árbol limpio; no modificar ramas protegidas.
- [x] 1.2 Registrar compileall, collect, Ruff y baseline de suite antes de código; documentar timeout/recuperación por particiones con límites explícitos.
- [ ] 1.3 Capturar galerías exploratorias actuales en `/tmp` y separar observaciones confirmadas de hipótesis.

## 2. Clean2D: ciclo bencénico manual
- [ ] 2.1 Reproducir un ciclo de seis átomos dibujado manualmente, compararlo con la ruta de plantilla y clasificar quick/publication/propose.
- [ ] 2.2 Registrar candidatos, fuentes, rechazos, IDs de selección, snapshot de grafo, longitudes, cruces, degeneración y undo.
- [ ] 2.3 Si se confirma deformación, añadir primero un test red/green que falle en el estado base y preserve enlaces/conectividad; corregir sólo la causa raíz.
- [ ] 2.4 Probar controles acíclicos, aromáticos explícitos, anillos simples, selecciones parciales, estilos estéreo y todas las invariantes.

## 3. Cuña, etiquetas y bifurcaciones
- [ ] 3.1 Extender la galería con F/Cl/Br, etiquetas de anchura variable y ramas de bifurcación.
- [ ] 3.2 Medir clearance etiqueta↔cuña y registrar la geometría antes del cambio en escala nominal/claro/oscuro.
- [ ] 3.3 Corregir trimming/encuentro únicamente si hay defecto confirmado; conservar API, orientación, recuento hash y sentido estéreo.
- [ ] 3.4 Añadir pruebas geométricas de no intersección, clearance acotado, continuidad de ramas y extremos/márgenes.

## 4. Canvas y anillos/Haworth
- [ ] 4.1 Reproducir qué límite impide crear enlaces junto al borde y rastrear eventos, coordenadas y `sceneRect`.
- [ ] 4.2 Habilitar expansión/navegación fuera del viewport visible sin clamp, mutación de grafo ni salto del viewport.
- [ ] 4.3 Revisar render, colocación, snapping y Undo de plantillas de anillo/Haworth; corregir sólo fallos reproducibles.
- [ ] 4.4 Añadir pruebas de interacción Qt con casos próximos al borde y ciclos/Haworth.

## 5. Validación y entrega
- [ ] 5.1 Producir una única lámina comparativa versionada y galería reproducible; dejar imágenes exploratorias en `/tmp`.
- [ ] 5.2 Ejecutar pruebas focales por fase, arquitectura, compileall, Ruff focal/global y OpenSpec estricto.
- [ ] 5.3 Ejecutar plan de ocho shards con límites por proceso; documentar cualquier timeout y reconciliar todos los node IDs sin declarar verde lo no verificado.
- [ ] 5.4 Actualizar `validation.md`, hacer commits pequeños y publicar sólo `feat/molecular-depiction-excellence`.
- [ ] 5.5 Informar SHA, commits, archivos, resultados, evidencia y estado `READY FOR MANUAL REVIEW` o `NOT READY`; no integrar/tag/release.
