# Tareas — delimitación del AI Molecular Structure Bridge

## Entregables de esta campaña de planificación

- [x] Confirmar `origin/main` vigente, baseline, limpieza del árbol y crear rama independiente.
- [x] Leer políticas de campañas, arquitectura y OpenSpecs de referencia.
- [x] Inspeccionar parser SMILES/ChemIO, MolGraph, ruta GUI de importación/inserción, Clean2D y resolver M16 name2structure.
- [x] Definir la arquitectura futura, interfaz de proveedor y contrato mínimo de solicitud/respuesta estructurada.
- [x] Definir validación ChemIO, estados/códigos estables, observabilidad, privacidad, límites y no-mutation.
- [x] Definir dependencias/prohibiciones y las pruebas deterministas offline requeridas.
- [x] Crear y validar este cambio OpenSpec de forma estricta.
- [x] Confirmar que no se cambió código, tests de producto, Clean2D ni `architecture/modules.yml`.

## Implementación futura Phase 1 — explícitamente NO ejecutada por este cambio

- [ ] Crear el módulo M23 propuesto en `src/chemuson/molecular_assistant/` (confirmar nombre/ruta), registrar ownership/API/dependencias en `architecture/modules.yml` y actualizar el contrato canónico de límites.
- [ ] Implementar tipos de solicitud, proveedor neutral, salida estructurada estricta `{ "smiles": "..." }` y resultados/razones estables.
- [ ] Implementar un único adaptador HTTP no streaming OpenAI-compatible Chat Completions, con endpoint/modelo explícitos, límites, timeout y secretos inyectados; no añadir soporte específico de llama.cpp/LM Studio.
- [ ] Validar mediante la ruta aislada existente `chemuson.chemio.rdkit_safe.smiles_to_molgraph_isolated`, sin fallback directo para datos no confiables; no crear parser ni importar GUI/Clean2D/ChemName.
- [ ] Añadir tests fake sin Internet/API para: solicitud vacía/excesiva; SMILES válido; SMILES inválido; JSON malformado; respuesta vacía o truncada; texto/fences/campos extra; timeout/cancelación explícita; error HTTP/red; respuesta oversized; códigos/estados estables; parsing determinista; equivalencia del grafo con importación ChemIO directa; distinción entre SMILES inválido y worker/parser no disponible; y ejecución sin credenciales.
- [ ] Añadir pruebas arquitectónicas de aislamiento de imports (sin GUI/Clean2D/ChemName en M23; sin M23/proveedores en M02) y ausencia de network en tests.
- [ ] Verificar que los errores producen resultados sin `MolGraph`; probar no-mutación de documento/canvas, selección, undo y dirty state con la integración gráfica correspondiente. La integración/atomicidad visual es gate del cambio futuro de UI, no de esta campaña.
- [ ] Ejecutar baseline y gates completos de arquitectura, suite, Ruff y OpenSpec en el cambio de implementación, investigar fallos nuevos y documentar fallos históricos.

## Cobertura mínima a conservar

Las pruebas futuras deben incluir expresamente los 12 casos solicitados: (1) fake válido, (2) SMILES inválido, (3) respuesta estructurada malformada, (4) timeout/red, (5) respuesta vacía, (6) texto adicional, (7) decodificación determinista, (8) ningún fallo muta la molécula existente, (9) sin imports GUI innecesarios, (10) Clean2D sin dependencia IA, (11) provider fake sin Internet y (12) construcción química idéntica a la importación SMILES ChemIO equivalente. Los casos (8) y la inserción normal sólo pueden cerrarse al implementar la integración gráfica posterior; su separación evita adelantar un cambio de UI a Phase 1.

## Fases posteriores — roadmap, no tareas de este cambio

1. Phase 2: comando/UI mínimo de dibujo a partir de descripción.
2. Phase 3: soporte de varios proveedores/modelos.
3. Phase 4: acciones moleculares estructuradas con validación/autorización.
4. Phase 5: evaluación Clean2D en OpenSpec separado; generar → validar → ejecutar Clean2D → comparar métricas. Sin dependencia M02→IA.

Permanecen fuera del plan Phase 1 los agentes autónomos, tool calling general, ejecución de código, edición directa del canvas/MolGraph, conversación persistente, memoria, RAG, búsqueda web, nomenclatura completa, mecanismos de reacción, modificación iterativa, optimización con feedback Clean2D, multimodalidad, integración ChemName, Campaign 5–9, rediseño UI, configuración completa de API keys y benchmarking.