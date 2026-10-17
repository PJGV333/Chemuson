# Tareas — CI y proyecciones

## 1. Baseline y OpenSpec
- [x] 1.1 Inspeccionar la especificación anterior, arquitectura, workflow y shard tool.
- [x] 1.2 Capturar baseline local y artefactos/logs de GitHub Actions antes de implementar.
- [x] 1.3 Validar esta propuesta estrictamente antes de tocar código.

## 2. CI con pesos por duración
- [x] 2.1 Añadir perfil versionado de tiempos reales por archivo y algoritmo LPT determinista con fallback.
- [x] 2.2 Instrumentar duración por test/archivo en artefactos aun ante timeout; preservar cobertura/terminación fail-closed.
- [x] 2.3 Actualizar workflow, contrato, tests, perfil y manifiesto exacto de node IDs; conservar timeout 300 s salvo evidencia justificada.
- [x] 2.4 Ejecutar ocho shards ponderados; la cobertura/JUnit/timings es exacta (2.358 IDs), con una falla basal de lifecycle Qt documentada en validation/AGENT_REPORT.

## 3. Viewport del Visual Lab
- [x] 3.1 Centrar canvas del laboratorio localmente sobre la geometría molecular tras show diferido y resize.
- [x] 3.2 Probar visibilidad real de átomos/enlaces/pixeles y sincronía Current/Working tras cambios de caso, zoom, rotación y resize.

## 4. Geometría y semántica de proyección
- [x] 4.1 Refinar Haworth canónica para piranosa/furanosa y posiciones por vértice, con validación de parámetros.
- [x] 4.2 Añadir descriptores tipados Haworth/Fischer/silla y adaptador de instancia seguro sin modificar MolGraph/.cmsn.
- [x] 4.3 Probar patrones β-D-glucopiranosa, inversión α/β y D/L, convención Fischer, relaciones axial/ecuatoriales/up/down y ausencia de stereo inventado.
- [ ] 4.4 Validar contactos con el renderer real a varias escalas y generar evidencia visual anterior/nueva a igual escala.

## 5. Validación y publicación
- [x] 5.1 Ejecutar focales, arquitectura, compileall, Ruff y OpenSpec estricto.
- [ ] 5.2 Ejecutar/verificar los ocho shards y CI completo en GitHub; corregir sin omitir tests.
- [ ] 5.3 Completar `validation.md`, commits pequeños, push a la misma rama; confirmar SHA remoto/local y worktree limpio. No hacer merge.
