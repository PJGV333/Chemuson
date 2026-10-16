# Tareas — excelencia geométrica de uniones

## 1. Rama, alcance y baseline
- [x] 1.1 Confirmar árbol limpio, fetch, SHA de prep y ramas protegidas sin modificarlas.
- [x] 1.2 Crear/publicar `feat/molecular-depiction-excellence` desde el SHA exacto de `origin/release/v0.3.0-beta.1-prep`.
- [x] 1.3 Capturar compileall, collect, Ruff, test focal y suite baseline por ocho shards; registrar F401 histórico sin corregir.
- [x] 1.4 Leer contratos GUI/estéreo previos y revisar `BondItem`, adaptadores de cuña, vecinos y estilo sin editar producción.

## 2. Diagnóstico y galería
- [x] 2.1 Añadir un generador reproducible de galería Qt con uniones simples/bold/múltiples, solid/hash y mixta, bifurcaciones 2–4, anillos aromáticos/fusionados, etiquetas/cargas/H, numeración, ángulos y grosores/escala diversos.
- [x] 2.2 Capturar la base en claro/oscuro, guardar una lámina comparativa reproducible y separar el saliente hash confirmado de la hipótesis sigma/pi; las capturas temporales individuales permanecen en `/tmp`.
- [x] 2.3 Medir perfiles sigma/pi y contornos hash en el encuentro doble/cuñas; no se confirmó un defecto sigma/pi visible, por lo que éstos no se cambiaron.
- [x] 2.4 Conservar los casos correctos como controles; no homogeneizar siluetas ni alterar proyección estereoquímica.

## 3. Geometría de uniones (condicionada a evidencia)
- [x] 3.1 Extraer el cálculo puro y acotado de segmentos hash a `src/chemuson/gui/bond_junctions.py`; la API de `BondItem` permanece igual.
- [x] 3.2 Evaluar perfiles sigma/pi; la hipótesis no quedó confirmada, así que no se añadió otro motor de perfiles ni se colapsaron dobles.
- [x] 3.3 No aplicar recortes/intersecciones no sustentados por geometría visible; mantener la construcción actual de solid wedge y sigma/pi.
- [x] 3.4 Centrar los hash marks en intervalos iguales y crecer el ancho linealmente en sus centros; la comparación visual del agente sugiere mejora en claro/oscuro, aceptación del propietario pendiente.
- [x] 3.5 Añadir pruebas geométricas para 30°, 60°, 90°, 120°, 150° y casi colineales en tres longitudes/grosores; la galería añade ramas, heteroátomos, cargas, H, numeración, anillos aromáticos y fusionados.
- [x] 3.6 Mantener sin cambio cualquier geometría sin estrategia confirmada y registrar sigma/pi como hipótesis pendiente.

## 4. Estilo Publication (gate explícito)
- [x] 4.1 Evaluar el gate: no se añade `ChemUSON Publication` porque no están validados todos los perfiles/exportaciones ni las preferencias/documentos.
- [x] 4.2 No modificar ChemDoodle-like, presets, preferencias ni serialización.
- [x] 4.3 Documentar la fase 4 como aplazada, no como terminada.

## 5. Validación y publicación
- [x] 5.1 Ejecutar pruebas focales (74), arquitectura/contrato CI (284), renderizado, doble enlace, exportación, templates/Haworth y aromáticos (74); la matriz cubre 30–179°, tres escalas y anchuras, y ambas apariencias.
- [x] 5.2 `compileall` completo y OpenSpec estricto pasan; Ruff global conserva sólo el F401 histórico registrado arriba. Se recolectaron 2267 IDs: siete shards terminaron en sus límites; el shard 4 agotó 300 s durante ejecución local simultánea. Se recuperó por familias acotadas: los 80 Clean2D IDs y 139 restantes (138 pass/1 skip) pasaron, y `test_ui_tool_rail.py` pasó 64/64 en aislamiento; cobertura global resultante 2247 pass, 20 skip, 0 fallos.
- [x] 5.3 Guardar [`docs/visual/bond-junctions-before-after.png`](../../../docs/visual/bond-junctions-before-after.png), lámina comparativa antes/después en claro/oscuro; mejora limitada y residuos/hipótesis documentados.
- [x] 5.4 Actualizar `validation.md`, crear commits descriptivos y publicar sólo `feat/molecular-depiction-excellence`; nunca integrar a prep/main/tag/release.
- [x] 5.5 Informar SHA, archivos, commits, pruebas/CI y `READY FOR MANUAL VISUAL REVIEW`, distinguiéndolo de aceptación visual del propietario.
