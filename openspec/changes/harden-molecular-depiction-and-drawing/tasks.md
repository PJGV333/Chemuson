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
- [x] 4.1 Reproducir por QTest que `_is_on_paper` rechazaba x=20 y x=-40 (margen 40); `sceneRect` no incluía coordenadas del grafo.
- [x] 4.2 Quitar el gate de página; expandir `sceneRect` por átomos + viewport actual y preservar el foco al Undo, sin clamp de coordenadas.
- [x] 4.3 Revisar Haworth α/β y front bold/plain; corregir la etiqueta OH de C3 que cruzaba dos enlaces del anillo y la etiqueta CH₂OH, preservando grafo y Undo.
- [x] 4.4 Añadir QTests de enlace en margen/borde/fuera del viewport, ciclo hexagonal en borde, expansión/Undo de sceneRect y colocación/Undo de Haworth.

## 5. Validación y entrega
- [x] 5.1 Versionar una única lámina before/after claro/oscuro con galería reproducible; dejar capturas exploratorias en `/tmp`.
- [x] 5.2 Ejecutar pruebas focales por fase, arquitectura, compileall, Ruff focal/global y OpenSpec estricto; documentar el hallazgo Ruff baseline conocido.
- [x] 5.3 Ejecutar plan de ocho shards con límites por proceso; shard 4 agotó 300 s y sus 284 IDs se verificaron por particiones. La ejecución oficial no se declara verde.
- [x] 5.4 Actualizar `validation.md` y guardar el trabajo en commits locales sólo en `feat/molecular-depiction-excellence`; no se publica sin autorización explícita.
- [x] 5.5 Informar SHA, commits, archivos, resultados, evidencia y estado `NOT READY`; no integrar/tag/release.
