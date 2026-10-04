# M23 — Asistente de estructuras moleculares

## Responsabilidad

M23 convierte una descripción textual en una propuesta estructurada `{"smiles":"..."}` y sólo publica un `MolGraph` después de que ChemIO acepte el SMILES mediante su worker aislado. El parser confirma aceptación sintáctica/química, no la correspondencia semántica entre la petición y la estructura propuesta.

## Límites

- Depende únicamente de M00 `core`, M01 `chemio` y la biblioteca estándar.
- No importa Clean2D (M02), ChemName (M04), GUI/controllers/canvas (M08–M13), name2structure (M16), composition root (M19) ni `tools`.
- M00, M01 y Clean2D no dependen de M23; Clean2D mantiene su pipeline determinista.
- El servicio no recibe documentos ni canvas, no persiste prompts/respuestas y no expone un grafo en resultados fallidos.
- El adaptador HTTP no tiene endpoint predeterminado: requiere base URL/modelo explícitos, usa Chat Completions no streaming, limita tiempo/tamaño y no sigue redirects.
- Las pruebas usan proveedores y transportes falsos; ninguna hace solicitudes de red.

## API pública

`MolecularAssistant`, `MolecularAssistantRequest`, `MolecularAssistantResult`, `MolecularAssistantStatus`, `MolecularStructureProvider`, `ProviderResponse`, `OpenAICompatibleConfig`, `OpenAICompatibleProvider`, `ProviderError`, `ProviderErrorCode` y `ProviderCancelled` se reexportan desde `chemuson.molecular_assistant`.
