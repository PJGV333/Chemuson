# Diseño: AI Molecular Structure Bridge

## Estado y límites

Este diseño describe Phase 1 / Foundation como trabajo futuro. El cambio activo actual sólo fija contratos; no implementa el servicio, transporte, adaptador, UI ni pruebas de runtime. El servicio no pertenece a Clean2D y no es un agente autónomo.

## Reconocimiento de ChemUSON

- **Modelo:** `src/chemuson/core/model.py` contiene `MolGraph` (M00), independiente de GUI.
- **Importador SMILES ordinario:** `chemuson.chemio.rdkit_io.smiles_to_molgraph` (M01) rechaza entrada vacía, intenta conversión a MOL mediante `rdkit_safe`/`_rdkit_worker.py` en un subprocess con timeout de 8 s y parsea el MOL con `molfile_to_molgraph`; si falla el worker, su comportamiento actual permite fallback a RDKit directo. Devuelve `MolGraph` o lanza una excepción.
- **Ruta aislada para salida no confiable:** `chemuson.chemio.rdkit_safe.smiles_to_molgraph_isolated` devuelve `(MolGraph | None, error | None)`, usa el worker existente (timeout de 8 s por defecto) y convierte el MOL mediante `molfile_to_molgraph`, sin fallback a RDKit directo en el proceso principal. Phase 1 debe usar esta ruta fail-closed para el SMILES generado y mapear sus fallos a estados/códigos estables. No se crea otro parser. Para una entrada aceptada, la prueba compara la identidad química resultante con `rdkit_io.smiles_to_molgraph`; ambos comparten el mismo worker/MOL parser en el camino normal. Si el worker no está disponible, el asistente devuelve un error de validación controlado en lugar de recurrir a la ruta directa menos aislada.
- **Depiction Clean2D:** `clean2d.imported_depiction.smiles_to_molgraph_best_depiction[_with_report]` selecciona geometría/candidatos usando M02 y no es validador autoritativo del SMILES. La GUI actual intenta esta ruta primero en `TemplateController._import_smiles_graph`, con fallback a ChemIO. El servicio IA no debe llamarla.
- **Inserción visual:** `TemplateController.on_import_smiles` entrega el grafo aceptado a `canvas._insert_molgraph`, que crea átomos/enlaces con comandos undoables; es un detalle de GUI/canvas y queda fuera de Phase 1. Una falla previa a la inserción debe producir cero cambios; la futura integración debe conservar undo/redo y probar la atomicidad de su commit.
- **Servicio existente distinto:** M16 `name2structure` resuelve nombres vía tabla estática y PubChem, valida por `rdkit_safe.smiles_to_molgraph_isolated` y retorna `NameToStructureResult`. No modela instrucciones en lenguaje natural ni proveedores sustituibles; no se convierte en el módulo IA.
- **Arquitectura:** M02 depende de M00/M01, M01 depende de M00 y ninguno necesita IA. M08/M10 son capas de orquestación de UI; M19 es sólo composition root de arranque, no un service locator ni un lugar para lógica de dominio.

## Módulo y dependencias propuestas

Reservar el siguiente ID disponible, M23, para el módulo de aplicación propuesto `molecular_assistant` en `src/chemuson/molecular_assistant/`. Mantener en ese módulo el caso de uso, protocolo de proveedor, adaptador OpenAI-compatible futuro, parser de respuesta y modelos de resultado: no crear una jerarquía de paquetes/proveedores ni un SDK/plugin framework para un solo proveedor.

Dependencias de M23 propuestas: M01 `chemio` para el worker/parser aislado existente, M00 `core` para exponer `MolGraph` y biblioteca estándar para JSON/contratos. M23 no importa M02, M04/ChemName, M08–M13/GUI, M16/name2structure ni M19. Los consumidores GUI futuros dependerán de M23; ni M00, M01 ni M02 dependerán de él. Actualizar el catálogo de módulos, la especificación de límites y tests arquitectónicos sólo en el cambio futuro que cree el paquete, evitando tocar `architecture/modules.yml` durante esta planificación.

```text
futura acción GUI / controller (fuera de Phase 1)
                    │ depende de
                    ▼
M23 Molecular Assistant ──► protocolo de proveedor
          │                          │
          │                          └─ futuro adaptador OpenAI Chat Completions
          ▼
respuesta JSON estricta {"smiles": "..."}
          │
          ▼
M01 chemio.smiles_to_molgraph ──► M00 MolGraph
          │
          └─ resultado tipado: estado + procedencia + validación
                                      │
                     GUI/canvas posterior: inserción estándar
                                      │
                           Clean2D opcional downstream
```

No se permite el flujo `Clean2D → IA`, ni que M23 importe/cargue la GUI o modifique un documento/canvas.

## Contratos de aplicación

### Solicitud y proveedor

- Solicitud mínima: un texto de usuario (`description`/prompt) no vacío y limitado en tamaño; no se incluye historial, memoria, selección del canvas, documento ni comandos de edición.
- Protocolo neutral al proveedor: operación única de generación que recibe la solicitud y el contrato/formato esperado, y devuelve contenido de respuesta como texto junto a `provider_id` y `model_id` opcional provenientes del transporte/configuración, nunca inferidos del JSON del modelo.
- Primer adaptador futuro: HTTP no streaming compatible con `POST /v1/chat/completions`, modelo configurado explícitamente y respuesta de contenido JSON. La interfaz de aplicación no depende de SDK, esquema de error OpenAI ni clases del transporte. No se da por hecho que todos los servidores implementen extensiones no comunes; verificar compatibilidad con un endpoint de referencia en la fase de implementación.
- No hay tool calling, retries de agente, streaming, funciones dinámicas, APIs de modelos ni petición directa desde GUI. Las credenciales, si una futura composición las inyecta, se tratan como secreto opaco, no se incluyen en prompts/logs/resultados; la UI/configuración persistente de claves queda fuera.

### Representación generada

La respuesta JSON de aplicación tiene exactamente un campo requerido: `smiles`, string no vacío. Ejemplo:

```json
{"smiles":"CCO"}
```

No se exige `name`: el usuario ya expresó qué solicita y la aplicación no debe confundir una etiqueta generada con una identidad química. No se exige `metadata`: proveedor/modelo son metadatos de transporte fiables y viven fuera del contenido del modelo. Campos adicionales se rechazan (no se ignoran silenciosamente ni se convierten en instrucciones). No aceptar Markdown, fences, prosa antes/después, texto truncado ni reparación automática del JSON.

### Validación y grafo

1. El adaptador extrae únicamente el campo de contenido definido por el protocolo de proveedor y aplica límites de bytes/tiempo antes de decodificarlo.
2. El decoder JSON acepta un objeto conforme exactamente al esquema anterior y límites de tamaño; cualquier respuesta incompleta, vacía, JSON inválido, forma/tipo incorrectos, campo desconocido o texto extra es `malformed_response`.
3. El caso de uso entrega la cadena `smiles` a `chemuson.chemio.rdkit_safe.smiles_to_molgraph_isolated`; no llama al parser de texto del proveedor, a M02 ni al canvas. Rechazo químico explícito del parser es `invalid_structure`; timeout/fallo del worker o de su infraestructura es `validation_error`, no un rechazo químico. Un SMILES aceptado sólo demuestra que el parser existente lo acepta: no demuestra que la estructura corresponda a «vancomicina» ni que el modelo sea científicamente correcto.
4. Sólo después de validar se publica un `MolGraph` en un resultado `success`. Los resultados que no son éxito no contienen grafo. La UI futura aplica ese grafo por sus mecanismos normales sólo tras recibir éxito, en una acción explícita y undoable.

### Resultado y diagnóstico estable

Cada resultado lleva `status` de vocabulario estable:

- `success`: parser aceptó SMILES y existe grafo;
- `invalid_request`: solicitud vacía o fuera de límites, rechazada antes de contactar al proveedor;
- `invalid_structure`: el parser ChemIO rechazó la estructura;
- `validation_error`: el worker/parser no pudo completar la validación por timeout o fallo de infraestructura;
- `provider_error`: timeout, red, rechazo HTTP/protocolo u otro error de proveedor;
- `malformed_response`: respuesta vacía/incompleta, JSON/schema inválidos, campos extra o contenido externo;
- `cancelled`: cancelación explícita antes de tener una estructura validada (si la implementación expone cancelación; nunca se presenta una respuesta parcial).

Además se devuelve `reason_code` estable cuando no hay éxito (por ejemplo `empty_prompt`, `request_too_large`, `timeout`, `network_error`, `http_error`, `invalid_json`, `unexpected_fields`, `missing_smiles`, `empty_smiles`, `response_too_large`, `invalid_smiles`, `parser_timeout`, `parser_unavailable`, `cancelled`), separado del diagnóstico técnico no estable. `cancelled` se emite sólo ante cancelación explícita; Phase 1 no añade un comando interactivo de cancelación. No traducir detalle crudo de excepciones al contrato ni usarlo como código. Éxito incluye `validation_passed=true`; una solicitud rechazada, error de proveedor/respuesta o incapacidad de validar antes de una decisión química usa `null`; rechazo químico usa `false`.

El resultado conserva `provider_id`, `model_id` opcional, el SMILES propuesto (si se pudo extraer), estado/validación y código estable. Un SMILES rechazado se conserva para diagnóstico sólo en el objeto de resultado durante la operación; no se escribe a documentos, historial, archivos de log ni telemetría persistente por defecto. Nunca guardar prompt, texto bruto, credenciales, cabeceras ni cuerpos HTTP.

## Seguridad, recursos y no-mutation

- Tratar solicitud, JSON, SMILES, errores y metadatos del modelo como datos no confiables. Parser de esquema estricto, sin coerción de tipos, extracción heurística, autocorrección ni ejecución de código.
- Contrato inicial recomendado: solicitud ≤ 8 KiB, cuerpo HTTP de respuesta ≤ 1 MiB, contenido JSON del modelo ≤ 16 KiB y SMILES ≤ 16 KiB; límites centralizados, comprobados antes de parsear y sin truncar. Ajustes requieren evidencia y pruebas, no ampliación ilimitada.
- Toda petición de proveedor tiene timeout finito y configurable; proponer 60 s como deadline inicial configurable de aplicación, con timeout de conexión acotado por el adaptador. La validación IA usa la ruta aislada existente de ChemIO, con su timeout de worker (8 s por defecto); si el worker falla, no se aplica fallback directo. Timeout, respuesta parcial, tamaño excesivo y error de red son resultados controlados, no excepciones que llegan a UI ni reintentos implícitos.
- Usar únicamente un endpoint configurado explícitamente; no derivar URL desde el prompt, seguir instrucciones del contenido o permitir que el modelo elija herramientas. No establecer conexión por defecto a servicio externo. No incluir secretos en errores/logs.
- El servicio no recibe documento/canvas y no importa GUI. Un fallo nunca produce `MolGraph` utilizable ni dispara callback de inserción; la molécula, selección, undo stack y dirty state actuales se mantienen intactos. La validación ocurre antes de commit. La atomicidad ante un error de aplicación visual se probará al planear la UI posterior.

## Observabilidad

Exponer al llamador campos estructurados: `provider_id`, `model_id` cuando esté disponible, `proposed_smiles` si fue extraído, `status`, `validation_passed` y `reason_code`. Mensajes/detalles técnicos son diagnósticos no estables y no deben sustituir códigos. Por defecto no persistir prompts, respuestas JSON completas, datos de usuario, credenciales ni SMILES en logs; el llamador decide si muestra el SMILES de la operación y aplica cualquier consentimiento/política de retención en otra fase.

## Pruebas previstas para Phase 1 y el límite UI

Usar fakes deterministas para el proveedor y el transporte; ninguna suite requiere API real, modelo local ni Internet. Probar respuesta válida (y equivalencia química con el parser ChemIO directo), SMILES inválido, JSON malformado, respuesta vacía/truncada, texto extra/fences/campos desconocidos, timeouts/error de red, límites, estados/códigos estables, parsing repetible, importación aislada sin GUI y ausencia de dependencias IA en Clean2D. Ver también el checklist de cobertura en `tasks.md`.

La prueba de que un fallo no modifica el canvas debe ejecutarse con el primer adaptador gráfico futuro, fuera de Phase 1: fake del provider + snapshot del grafo/selección/undo/dirty-state antes y después del error, verificando que no se llama a inserción. Phase 1 prueba que el servicio carece de acceso/capacidad para hacer esa mutación.

## Riesgos y decisiones pospuestas

- El importador ordinario `smiles_to_molgraph` hace fallback a RDKit en proceso si falla el worker; la frontera IA elige intencionalmente el helper aislado ya existente y falla de forma controlada en esa situación. La equivalencia química se prueba cuando ambas rutas aceptan la entrada; no se introduce otro parser ni se altera el importador ordinario.
- La compatibilidad real de `response_format`/JSON mode varía entre servidores “OpenAI-compatible”; el contrato valida JSON estrictamente aunque el servidor no garantice modo JSON. La primera implementación probará al menos un endpoint local compatible; llama-server y LM Studio son ejemplos de servidores, no dependencias de producción.
- La inserción SMILES actual selecciona candidatos M02 antes de ChemIO. El módulo IA no debe llamar esa fachada para validación; futura integración decidirá si representa el grafo aceptado tal cual o invoca Clean2D separadamente como operación gráfica explícita.
- El Phase 1 no afirma exactitud química semántica, rendimiento de modelos ni mejora de Clean2D.

## Roadmap no vinculante

1. **Phase 1 / Foundation:** servicio, protocolo, un adaptador compatible OpenAI Chat Completions, salida JSON estricta con SMILES, validación ChemIO, diagnóstico y tests fake/offline.
2. **Phase 2:** comando/UI mínimo «dibujar molécula a partir de descripción» e inserción segura/undoable.
3. **Phase 3:** más proveedores/modelos tras contratos y compatibilidad medida.
4. **Phase 4:** acciones moleculares estructuradas con autorización y validaciones adicionales.
5. **Phase 5:** campaña separada de evaluación Clean2D (estructura validada → Clean2D → métricas before/after), con su propio OpenSpec, datos y gates. Nunca una dependencia M02→IA.
