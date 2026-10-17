# Tareas — Preset AB3 y proyecciones especiales

## 1. Contrato y baseline
- [x] 1.1 Leer el OpenSpec del laboratorio, sus contratos, arquitectura M08 y constructores/fixtures actuales.
- [x] 1.2 Registrar baseline completo antes de modificar código o documentación.
- [x] 1.3 Validar esta propuesta OpenSpec estrictamente antes de implementar.

## 2. Preset y corpus editorial
- [ ] 2.1 Implementar el preset nombrado `ChemUSON Editorial AB3` con los parámetros exactos y sin tocar defaults globales.
- [ ] 2.2 Permitir cargarlo desde el laboratorio y CLI/exportaciones reproducibles.
- [ ] 2.3 Marcar elegibilidad y rol por fixture; excluir B4 de la selección editorial sin quitarlo de render.
- [ ] 2.4 Aclarar elegibilidad/exclusión en exportes, resumen y contacto de galería.
- [ ] 2.5 Añadir pruebas de valores exactos, Current inmutable, carga y exclusión B4.

## 3. Geometría de proyecciones
- [ ] 3.1 Añadir tipos geométricos puros compartidos y registrar el subpaquete bajo M08.
- [ ] 3.2 Reconstruir Haworth desde parámetros/reglas explícitas y conservar adaptador heredado.
- [ ] 3.3 Añadir parámetros/constructores extensibles para Fischer y silla.
- [ ] 3.4 Probar determinismo, topología, frontales, stereo no inventado, colisiones y límites de imports.

## 4. Evidencia y cierre
- [ ] 4.1 Generar comparación visual anterior/nueva de Haworth y exportación editorial AB3 reproducible.
- [ ] 4.2 Ejecutar pruebas focales, arquitectura, compileall, Ruff focal y shards completos acotados.
- [ ] 4.3 Validar OpenSpec, diff y reconciliar el fallo heredado de baseline sin ocultarlo.
- [ ] 4.4 Completar `validation.md` con archivos, decisiones, evidencia y límites.
