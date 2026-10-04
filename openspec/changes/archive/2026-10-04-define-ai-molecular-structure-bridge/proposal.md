# Propuesta: frontera para asistencia molecular mediante IA

## Propósito

Definir una frontera de aplicación pequeña, segura y desacoplada para que ChemUSON pueda solicitar a un proveedor compatible con OpenAI una propuesta de estructura molecular en SMILES, validarla con el parser ChemIO existente (priorizando su ruta aislada para salida no confiable) y devolver un `MolGraph` aceptado por ChemUSON. El modelo propone; el parser y el dominio de ChemUSON deciden. Esta capacidad constituye una campaña independiente y no forma parte de Clean2D.

La campaña implementa la Phase 1 / Foundation para interpretar una petición como «dibuja vancomicina» detrás de una frontera provider-neutral. No incorpora interfaz, no pretende demostrar la exactitud semántica de la propuesta y no pertenece a Clean2D.

## Alcance de este cambio

- Un servicio de aplicación molecular y un contrato de proveedor sustituible.
- Un adaptador HTTP no streaming para Chat Completions compatible con OpenAI, con endpoint/modelo explícitos, límites, timeout finito y credenciales opcionales inyectadas.
- Un contrato mínimo cuya única salida molecular requerida es `{"smiles":"..."}`, decodificada sin extracción ni reparación.
- Validación del SMILES mediante la ruta aislada existente de ChemIO y resultados/códigos estables.
- Límites de arquitectura, privacidad, observabilidad y pruebas deterministas offline con providers/transports falsos.

## Resultado de Phase 1

M23 (`src/chemuson/molecular_assistant/`) y sus pruebas están implementados y registrados en `architecture/modules.yml`. Se actualizaron los deltas OpenSpec de catálogo y límites para conservar la dirección de dependencias. No se añade una dependencia externa ni se establece una conexión de red durante importación o por defecto.

## Fuera de alcance

Quedan fuera la conexión a un modelo/API real durante tests o por defecto, configuración/UI de API keys, UI/command de usuario, llama.cpp/llama-server, integración específica de LM Studio, varios proveedores, agentes, tool calling, ejecución de código, modificación directa del canvas o de `MolGraph`, memoria, conversación persistente, RAG, búsqueda web, nomenclatura química completa, generación de reacciones, edición molecular iterativa, feedback/optimización Clean2D, multimodalidad, integración con ChemName, benchmarking, cambios grandes de UI y Campaign 5–9 de Clean2D. Tampoco se modifican algoritmos, scoring, candidatos, políticas o corpus/regresiones de Clean2D.

## Módulos y dirección de dependencias

- **M23 `molecular_assistant` (`src/chemuson/molecular_assistant/`):** posee la solicitud, el contrato/adaptador de proveedor, la decodificación estricta y el resultado trazable. No importa GUI ni Clean2D.
- **M01 `chemio`:** conserva el parsing/validación de SMILES existente y no adquiere dependencia de IA.
- **M00 `core`:** aporta el tipo de salida `MolGraph`.
- **M08/M10/M09 (`gui`, controllers, canvas):** consumidores futuros, fuera de Phase 1.
- **M02 `clean2d`:** permanece downstream/determinista y no importa ni contacta IA.

`architecture/modules.yml`, los tests de límites y los deltas OpenSpec registran M23 con dependencias únicamente M00/M01, y prohíben las dependencias inversas desde M00/M01/M02.

## Suposición arquitectónica comprobada

El repositorio ya expone `chemuson.chemio.rdkit_io.smiles_to_molgraph`, que intenta el worker RDKit aislado y tiene un fallback existente a RDKit directo. La UI de importación SMILES actual no es sólo ese parser: `TemplateController._import_smiles_graph` intenta primero selección de depiction de M02 (`smiles_to_molgraph_best_depiction_with_report`) y después usa el importador ChemIO; la inserción visual llama `canvas._insert_molgraph`. El servicio IA valida con ChemIO y no llama a la ruta de depiction Clean2D. La integración gráfica posterior deberá decidir explícitamente cómo reutilizar la inserción normal y si solicita Clean2D de manera opcional.

También existe M16 `name2structure`, un resolver de nombres comunes/PubChem con conectores propios, no un asistente lingüístico ni una interfaz de proveedor LLM. No se reutiliza ni se amplía en esta campaña.

## Criterio de éxito

Phase 1 ofrece un servicio reemplazable y comprobable offline, sólo publica un `MolGraph` validado, clasifica fallos con razones estables y mantiene aislados Clean2D y la GUI. Las pruebas de contrato, el catálogo arquitectónico y la validación estricta OpenSpec verifican esos límites.