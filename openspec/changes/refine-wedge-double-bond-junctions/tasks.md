# Tareas — refinamiento de wedges y dobles enlaces

## 1. Rama, baseline y OpenSpec
- [x] 1.1 Confirmar rama/HEAD/upstream limpios y conservar checkpoint local antes de producción.
- [x] 1.2 Capturar compileall, collect, Ruff selectivo y baseline pytest con timeout; registrar hallazgos existentes.
- [x] 1.3 Guardar render baseline actual en `/tmp/chemuson-depiction-refinement-baseline/before-existing/`.
- [x] 1.4 Leer contratos de unión/label existentes y registrar hipótesis de causa raíz sin parchear.

## 2. Reproducción y evidencia visual
- [ ] 2.1 Ampliar galería con hash↔doble y su línea pi próxima, usando la ruta de producción.
- [ ] 2.2 Añadir O/N/halógenos cortos junto a dobles horizontales, verticales y diagonales; medir centro óptico vs átomo y trazos.
- [ ] 2.3 Capturar wedge sólido/hash con etiquetas y medir gap desde la silueta real del trazo, incluidas caps.
- [ ] 2.4 Añadir Haworth pentagonal y hexagonal con aristas bold y wedges sólidos en claro/oscuro.
- [ ] 2.5 Confirmar cada defecto con comparación base y guardar crops/intermedios bajo `/tmp`.

## 3. Correcciones geométricas condicionadas a evidencia
- [ ] 3.1 Corregir intersección/trim hash↔doble sin confundir sigma con pi ni alterar hash count/orientación.
- [ ] 3.2 Mejorar el centrado óptico de etiquetas en dobles sólo en los contextos demostrados.
- [ ] 3.3 Aplicar política de clearance consistente con ancho/cap del trazo, gap mínimo y máximo medidos.
- [ ] 3.4 Refinar wedge sólido↔bold en Haworth si tests/galería confirman bulto; mantener cuña limitada y no autointersectada.
- [ ] 3.5 Conservar identidad del grafo, coordenadas, estéreo, semántica bold y los casos aceptados de la campaña previa.

## 4. Regresión y entrega
- [ ] 4.1 Añadir tests para hash↔doble, clearance wedge-label y wedge sólido↔bold pentágono/hexágono.
- [ ] 4.2 Ejecutar controles de doble normal, wedges sin label, heteroátomos, anillo sencillo, Haworth previo y galería anterior.
- [ ] 4.3 Actualizar la manifest CI sólo si cambian node IDs; revisar las adiciones antes de regenerar.
- [ ] 4.4 Ejecutar pruebas focales, arquitectura, compileall, Ruff focal/global, OpenSpec estricto y shards con límites; reconciliar timeouts honestamente.
- [ ] 4.5 Versionar before/after claro/oscuro, actualizar validation, realizar commits pequeños y publicar únicamente esta rama experimental.
- [ ] 4.6 Reportar SHA/commits, rutas de evidencia, resultado por prioridad y `READY FOR MANUAL REVIEW` o `NOT READY`; sin merge/tag/release.
