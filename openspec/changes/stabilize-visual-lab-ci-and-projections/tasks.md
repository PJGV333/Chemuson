# Tareas — CI y proyecciones

## 1. Baseline y OpenSpec
- [x] 1.1 Inspeccionar la especificación anterior, arquitectura, workflow y shard tool.
- [x] 1.2 Capturar baseline local y artefactos/logs de GitHub Actions antes de implementar.
- [x] 1.3 Validar esta propuesta estrictamente antes de tocar código.

## 2. CI con pesos por duración
- [ ] 2.1 Añadir perfil versionado de tiempos reales por archivo y algoritmo LPT determinista con fallback.
- [ ] 2.2 Instrumentar duración por test/archivo en artefactos aun ante timeout; preservar cobertura/terminación fail-closed.
- [ ] 2.3 Actualizar workflow, contrato, tests, perfil y manifiesto exacto de node IDs; conservar timeout 300 s salvo evidencia justificada.
- [ ] 2.4 Simular y ejecutar ocho shards ponderados, validar cobertura exacta y documentar causa raíz.

## 3. Viewport del Visual Lab
- [ ] 3.1 Centrar canvas del laboratorio localmente sobre la geometría molecular tras show diferido y resize.
- [ ] 3.2 Probar visibilidad real de átomos/enlaces/pixeles y sincronía Current/Working tras cambios de caso, zoom, rotación y resize.

## 4. Geometría y semántica de proyección
- [ ] 4.1 Refinar Haworth canónica para piranosa/furanosa y posiciones por vértice, con validación de parámetros.
- [ ] 4.2 Añadir descriptores tipados Haworth/Fischer/silla y adaptador de instancia seguro sin modificar MolGraph/.cmsn.
- [ ] 4.3 Probar patrones β-D-glucopiranosa, inversión α/β y D/L, convención Fischer, relaciones axial/ecuatoriales/up/down y ausencia de stereo inventado.
- [ ] 4.4 Validar contactos con el renderer real a varias escalas y generar evidencia visual anterior/nueva a igual escala.

## 5. Validación y publicación
- [ ] 5.1 Ejecutar focales, arquitectura, compileall, Ruff y OpenSpec estricto.
- [ ] 5.2 Ejecutar/verificar los ocho shards y CI completo en GitHub; corregir sin omitir tests.
- [ ] 5.3 Completar `validation.md`, commits pequeños, push a la misma rama; confirmar SHA remoto/local y worktree limpio. No hacer merge.
