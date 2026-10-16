# Tareas — robustez de depiction y dibujo

## 1. Rama y baseline
- [x] 1.1 Confirmar rama experimental/upstream, SHA inicial y árbol limpio; no modificar ramas protegidas.
- [x] 1.2 Registrar compileall, collect, Ruff y baseline de suite antes de código; documentar timeout/recuperación por particiones con límites explícitos.
- [x] 1.3 Capturar la galería actual en `/tmp/chemuson-depiction-campaign-baseline/current-gallery/`; registrar mediciones y distinguir observaciones reproducidas de hipótesis.

## 2. Clean2D: ciclo bencénico manual
- [ ] 2.1 Reproducir un ciclo de seis átomos dibujado manualmente, compararlo con la ruta de plantilla y clasificar quick/publication/propose.
- [x] 2.2 Registrar fuentes/rechazos, selección, snapshot de grafo, longitudes, cruces, degeneración y Undo en sondas del motor y GUI; artefactos temporales en `/tmp/chemuson-depiction-campaign-baseline/`.
- [ ] 2.3 Si se confirma deformación, añadir primero un test red/green que falle en el estado base y preserve enlaces/conectividad; corregir sólo la causa raíz.
- [ ] 2.4 Probar controles acíclicos, aromáticos explícitos, anillos simples, selecciones parciales, estilos estéreo y todas las invariantes.

## 3. Cuña, etiquetas y bifurcaciones
- [x] 3.1 Extender la galería con F/Cl/Br/N⁺/CH₂OH, punta y base, cuña/hash, estilos sigma/pi y bifurcaciones.
- [x] 3.2 Medir clearance nominal antes/después en ambos temas; artefactos de comparación temporales en `/tmp`.
- [x] 3.3 Corregir el exceso de trimming confirmado por medición geométrica; conservar API, orientación, recuento/hash marks y sentido estéreo.
- [x] 3.4 Añadir matriz geométrica de no intersección y clearance acotado para F/Cl/Br/N⁺/CH₂OH, 4 orientaciones, punta/base, simple/doble/triple/bold/cuñas.

## 4. Canvas y anillos/Haworth
- [ ] 4.1 Reproducir qué límite impide crear enlaces junto al borde y rastrear eventos, coordenadas y `sceneRect`.
- [ ] 4.2 Habilitar expansión/navegación fuera del viewport visible sin clamp, mutación de grafo ni salto del viewport.
- [ ] 4.3 Revisar render, colocación, snapping y Undo de plantillas de anillo/Haworth; corregir sólo fallos reproducibles.
- [ ] 4.4 Añadir pruebas de interacción Qt con casos próximos al borde y ciclos/Haworth.

## 5. Validación y entrega
- [x] 5.1 Versionar una única lámina before/after claro/oscuro con galería reproducible; dejar capturas exploratorias en `/tmp`.
- [ ] 5.2 Ejecutar pruebas focales por fase, arquitectura, compileall, Ruff focal/global y OpenSpec estricto.
- [ ] 5.3 Ejecutar plan de ocho shards con límites por proceso; documentar cualquier timeout y reconciliar todos los node IDs sin declarar verde lo no verificado.
- [ ] 5.4 Actualizar `validation.md`, hacer commits pequeños y publicar sólo `feat/molecular-depiction-excellence`.
- [ ] 5.5 Informar SHA, commits, archivos, resultados, evidencia y estado `READY FOR MANUAL REVIEW` o `NOT READY`; no integrar/tag/release.
