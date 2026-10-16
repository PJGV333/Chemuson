# Tareas — excelencia geométrica de uniones

## 1. Rama, alcance y baseline
- [x] 1.1 Confirmar árbol limpio, fetch, SHA de prep y ramas protegidas sin modificarlas.
- [x] 1.2 Crear/publicar `feat/molecular-depiction-excellence` desde el SHA exacto de `origin/release/v0.3.0-beta.1-prep`.
- [x] 1.3 Capturar compileall, collect, Ruff, test focal y suite baseline por ocho shards; registrar F401 histórico sin corregir.
- [x] 1.4 Leer contratos GUI/estéreo previos y revisar `BondItem`, adaptadores de cuña, vecinos y estilo sin editar producción.

## 2. Diagnóstico y galería
- [ ] 2.1 Añadir un generador reproducible de una galería Qt vectorial con uniones simples/bold/múltiples, solid/hash y mixta, bifurcaciones 2–4, anillos aromáticos/fusionados, etiquetas/cargas/H, ángulos y grosores/escala diversos.
- [ ] 2.2 Capturar la base en temas claro y oscuro y documentar defectos confirmados frente a hipótesis; no versionar capturas temporales.
- [ ] 2.3 Medir perfiles sigma/pi y contornos rayados en el encuentro doble/cuñas; confirmar si su tratamiento actual causa un defecto visible antes de editarlo.
- [ ] 2.4 Conservar los casos correctos como controles; no homogeneizar siluetas ni alterar proyección estereoquímica.

## 3. Geometría de uniones (condicionada a evidencia)
- [ ] 3.1 Si una insuficiencia geométrica queda confirmada, extraer a `src/chemuson/gui/bond_junctions.py` únicamente cálculos puros necesarios, manteniendo las API existentes de `BondItem`.
- [ ] 3.2 Representar dobles como perfiles sigma/pi independientes, con offsets, grosores y trims efectivos; no modelarlos como una sola línea ancha.
- [ ] 3.3 Aplicar recortes/intersecciones acotados, derivados de contornos/anchuras reales; no añadir listas nuevas de compensaciones arbitrarias.
- [ ] 3.4 Mejorar cuñas rayadas sólo si la distribución/contacto mejora visualmente en todo el rango de escalas y uniones.
- [ ] 3.5 Añadir casos de ángulos 30°, 60°, 90°, 120°, 150° y casi colineales, escalas/grosores, bifurcaciones, etiquetas/cargas/H, anillos y dobles aromáticos.
- [ ] 3.6 Si no existe estrategia segura, conservar el último estado válido, registrar el defecto/hipótesis pendiente y no cambiar geometría de producción.

## 4. Estilo Publication (gate explícito)
- [ ] 4.1 Añadir `ChemUSON Publication` como estilo opcional sólo si las fases 1–3 quedan validadas con imagen antes/después, suite relevante y contrato de preferencias/documentos.
- [ ] 4.2 Mantener ChemDoodle-like y todas las preferencias/documentos existentes sin cambio silencioso; comprobar retorno al estilo actual.
- [ ] 4.3 Si el gate visual o de preferencias no se cumple, documentar la fase 4 como aplazada, no como terminada.

## 5. Validación y publicación
- [ ] 5.1 Ejecutar pruebas focales, renderizado GUI, arquitectura y regresiones Haworth/aromáticos/complejas.
- [ ] 5.2 Ejecutar `compileall`, Ruff focalizado (sin perseguir F401 histórico), OpenSpec estricto, `git diff --check` y pytest completo por shards con límites externos.
- [ ] 5.3 Guardar una lámina comparativa antes/después, claro/oscuro, más una matriz concisa de mejoras confirmadas y residuos.
- [ ] 5.4 Actualizar `validation.md`, commits pequeños descriptivos y publicar sólo `feat/molecular-depiction-excellence`; nunca integrar a prep/main/tag/release.
- [ ] 5.5 Informar SHA, archivos, commits, pruebas/CI y `READY FOR MANUAL VISUAL REVIEW` o `NOT READY`, distinguiéndolo de aceptación visual del propietario.
