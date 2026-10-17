# Tareas — Preset AB3 y proyecciones especiales

## 1. Contrato y baseline
- [x] 1.1 Leer el OpenSpec del laboratorio, sus contratos, arquitectura M08 y constructores/fixtures actuales.
- [x] 1.2 Registrar baseline completo antes de modificar código o documentación.
- [x] 1.3 Validar esta propuesta OpenSpec estrictamente antes de implementar.

## 2. Preset y corpus editorial
- [x] 2.1 Implementar el preset nombrado `ChemUSON Editorial AB3` con los parámetros exactos y sin tocar defaults globales.
- [x] 2.2 Permitir cargarlo desde el laboratorio y CLI/exportaciones reproducibles.
- [x] 2.3 Marcar elegibilidad y rol por fixture; incluir sólo A1, A2, A3, B1, B2 y B3; excluir A4 del corpus aprobado y mantener B4 renderizable como fixture no editorial.
- [x] 2.4 Aclarar elegibilidad/exclusión en exportes, resumen y contacto de galería.
- [x] 2.5 Añadir pruebas de valores exactos, Current inmutable, carga y exclusiones A4/B4.

## 3. Geometría de proyecciones
- [x] 3.1 Añadir tipos geométricos puros compartidos y registrar el subpaquete bajo M08.
- [x] 3.2 Reconstruir Haworth desde parámetros/reglas explícitas y conservar adaptador heredado.
- [x] 3.3 Añadir parámetros/constructores extensibles para Fischer y silla.
- [x] 3.4 Probar determinismo, topología, frontales, stereo no inventado, colisiones y límites de imports.

## 4. Evidencia y cierre
- [x] 4.1 Generar comparación visual anterior/nueva de Haworth y exportación editorial AB3 reproducible.
- [x] 4.2 Ejecutar pruebas focales, arquitectura, compileall, Ruff focal y shards completos acotados.
- [x] 4.3 Validar OpenSpec, diff y reconciliar el fallo heredado de baseline sin ocultarlo.
- [x] 4.4 Completar `validation.md` con archivos, decisiones, evidencia y límites.
