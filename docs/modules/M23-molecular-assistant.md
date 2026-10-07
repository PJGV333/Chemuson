# M23 — Asistente de estructuras moleculares

## Responsabilidad

M23 convierte una descripción textual en una propuesta estructurada `{"smiles":"..."}` y sólo publica un `MolGraph` después de que ChemIO acepte el SMILES mediante su worker aislado. El parser confirma aceptación sintáctica/química, no la correspondencia semántica entre la petición y la estructura propuesta.

## Límites

- Depende únicamente de M00 `core`, M01 `chemio` y la biblioteca estándar.
- No importa Clean2D (M02), ChemName (M04), GUI/controllers/canvas (M08–M13), name2structure (M16), composition root (M19) ni `tools`.
- M10 `gui.controllers` puede consumir M23 como adaptador asíncrono para la UI; M23 no importa ni depende de GUI/controllers.
- M00, M01 y Clean2D no dependen de M23; Clean2D mantiene su pipeline determinista.
- El servicio no recibe documentos ni canvas, no persiste prompts/respuestas y no expone un grafo en resultados fallidos.
- El adaptador HTTP requiere base URL/modelo explícitos, usa Chat Completions no streaming, limita tiempo/tamaño y no sigue redirects. Los perfiles OpenAI, LM Studio y llama.cpp sólo aportan defaults editables; no añaden protocolos ni afirman compatibilidad live.
- Los IDs de modelo son aportados por el usuario porque dependen del endpoint cargado; no hay catálogo ni descubrimiento de modelos.
- Las pruebas usan proveedores y transportes falsos; ninguna hace solicitudes de red.

## API pública

`MolecularAssistant`, `MolecularAssistantRequest`, `MolecularAssistantResult`, `MolecularAssistantStatus`, `MolecularStructureProvider`, `ProviderResponse`, `OpenAICompatibleConfig`, `OpenAICompatibleProfile`, `OPENAI_COMPATIBLE_PROFILES`, `get_openai_compatible_profile`, `OpenAICompatibleProvider`, `ProviderError`, `ProviderErrorCode` y `ProviderCancelled` se reexportan desde `chemuson.molecular_assistant`.

## Limitaciones conocidas — integrable, pero no infalible

### Limitaciones de modelos y endpoints

- Un LLM puede producir un SMILES químicamente válido que represente una molécula distinta de la solicitada.
- El LLM puede responder con JSON inválido o no respetar el esquema textual estricto.
- El endpoint NInfer usado en las pruebas actuales no soporta `response_format={"type":"json_object"}`; ChemUSON usa el fallback textual estricto. Es una limitación de compatibilidad del endpoint, no evidencia de un fallo químico de ChemUSON.
- Los modelos razonadores pueden consumir `max_tokens` en reasoning y devolver `content` vacío. ChemUSON lo clasifica como `generation_exhausted`; no existe contenido que reparar y no se conserva `reasoning_content`.
- No se garantiza que IA solamente recupere moléculas complejas como tetrandrina, vancomicina o eritromicina.

### Límites de validación y resolución de ChemUSON

- ChemIO aislado valida que el SMILES produzca una estructura interpretable; por sí solo no demuestra que sea la molécula pedida.
- Identity/reference verification reduce el riesgo comparando identidad aislada con una referencia, pero el resultado puede permanecer `UNVERIFIED` si no hay referencia confiable disponible.
- El modo offline no consulta PubChem: permite referencias locales y la caché local existente. El acceso externo requiere permiso explícito y lo ejecuta ChemUSON, no el modelo.
- Whole-molecule transformations siguen siendo IA-only y no disponen necesariamente de una referencia determinista.
- Qwen no recibe navegador, herramientas web ni URL de PubChem; no realiza HTTP por su cuenta. Si el usuario permite lookup externo, ChemUSON envía a PubChem únicamente el nombre químico extraído.

### Capacidades todavía no ofrecidas

- No existe inicio de sesión GPT ni OAuth de cuenta.
- No se ofrece browsing del modelo (intencionalmente), un modelo local especializado en química ni fine-tuning químico.

### Deuda técnica conocida de ChemUSON

- Persiste una deuda Qt independiente: `_DescriptorWorker` puede provocar un fallo de teardown en ciertas ejecuciones agregadas de tests. El bug concreto del wrapper stale de `MolecularAssistantDialog` sí se corrigió; eso no demuestra que todos los SIGSEGV/teardown Qt estén resueltos.
- Sigue existiendo deuda global de OpenSpec: specs de capacidades no relacionadas conservan placeholders en su sección Purpose. La validación estricta de este cambio activo se informa por separado; no se modifican esos placeholders aquí.

La incompatibilidad detectada en el contrato de propiedades PUG REST de Name→Structure fue un bug del conector ChemUSON y se corrigió en el cierre de referencia; no debe atribuirse a una alucinación del modelo.
