# Tareas — AI Molecular Structure Bridge / Phase 1 Foundation

## Preparación OpenSpec completada antes de implementar

- [x] Confirmar `origin/main`, baseline documentado, árbol limpio y rama `ai/molecular-assistant-foundation`.
- [x] Leer políticas de campañas, arquitectura y OpenSpecs de referencia.
- [x] Inspeccionar parser SMILES/ChemIO, MolGraph, flujo GUI de importación/inserción, Clean2D y M16 `name2structure`.
- [x] Definir contrato, validación ChemIO, estados/códigos, observabilidad, privacidad, límites y no-mutation.
- [x] Validar los contratos OpenSpec estrictamente antes de implementar.

## Implementación Phase 1

- [x] Crear M23 `src/chemuson/molecular_assistant/`, registrar ownership/API/dependencias y fijar contratos de catálogo/límites.
- [x] Implementar request, provider-neutral protocol, salida JSON estricta `{ "smiles": "..." }` y resultados/códigos estables.
- [x] Implementar un adaptador HTTP no streaming OpenAI-compatible Chat Completions con endpoint/modelo explícitos, límites, timeout finito y secreto opcional inyectado; sin soporte específico de llama.cpp/LM Studio.
- [x] Validar el SMILES exclusivamente mediante `chemuson.chemio.rdkit_safe.smiles_to_molgraph_isolated`, sin fallback directo, parser nuevo, GUI, Clean2D ni ChemName.
- [x] Añadir tests offline/fake para request vacía/excesiva; estructura válida e inválida; JSON vacío/malformado/truncado/extra/fences/campos adicionales; límites; timeout/cancelación; errores HTTP/red; parsing determinista; equivalencia ChemIO cuando RDKit worker está disponible; errores exactos `invalid_input`/`timeout`/`rdkit_unavailable` y fallback genérico sin inspeccionar diagnósticos; y provider sin credenciales/red.
- [x] Añadir tests arquitectónicos para imports M23, aislamiento de GUI/Clean2D/ChemName, dependencia inversa prohibida desde M00/M01/M02, import aislado y ausencia de contexto documento/canvas en el request.
- [x] Verificar que todo resultado fallido carece de `MolGraph`; M23 no recibe documento/canvas y no puede efectuar inserción.
- [ ] Probar snapshot de grafo/selección/undo/dirty-state con integración visual; gate de una campaña UI posterior, fuera de Phase 1.
- [x] Ejecutar y registrar los gates finales: arquitectura, suite completa, Ruff, OpenSpec estricto, compileall, colección y `git diff --check`; documentar el único fallo histórico y el F401 preexistente.

## Cobertura y fases posteriores

La cobertura Phase 1 conserva los casos fake válido, SMILES inválido, respuesta malformada/vacía/texto extra, timeout/red, decodificación determinista, imports GUI aislados, Clean2D sin dependencia IA, fake sin Internet y equivalencia química con ChemIO cuando el worker está disponible. La no-mutación observable del canvas y la inserción normal se verificarán al implementar Phase 2.

1. **Phase 1 — Foundation:** completada por esta campaña.
2. **Phase 2 — UI/comando mínimo:** descripción → estructura validada → inserción normal y undoable.
3. **Phase 3 — evaluación Clean2D:** OpenSpec separado; la orquestación puede invocar por separado M23/M02, nunca crear M02 → M23.
4. **Phase 4 — más proveedores/modelos:** sólo tras contratos y compatibilidad medida.
5. **Phase 5 — edición molecular estructurada:** autorización y validaciones adicionales.

Siguen fuera de Phase 1 los agentes, tool calling, ejecución de código, edición directa del canvas/MolGraph, conversación/memoria/RAG, búsqueda web, nomenclatura completa, reacciones, edición iterativa, feedback Clean2D, multimodalidad, ChemName, Campaign 5–9, rediseño UI, configuración persistente de API keys y benchmarking.
