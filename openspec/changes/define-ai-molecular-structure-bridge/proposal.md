# Propuesta: frontera para asistencia molecular mediante IA

## Propósito

Definir una frontera de aplicación pequeña, segura y desacoplada para que ChemUSON pueda solicitar a un proveedor compatible con OpenAI una propuesta de estructura molecular en SMILES, validarla con el parser ChemIO existente (priorizando su ruta aislada para salida no confiable) y devolver un `MolGraph` aceptado por ChemUSON. El modelo propone; el parser y el dominio de ChemUSON deciden. Esta capacidad constituye una campaña independiente y no forma parte de Clean2D.

La campaña responde a una necesidad futura concreta —interpretar «dibuja vancomicina»— sin implementar todavía interfaz, transporte HTTP ni llamadas a modelos.

## Alcance de este cambio

Este OpenSpec delimita el Phase 1 / Foundation que se implementará en un cambio posterior:

- un servicio de aplicación molecular y un contrato de proveedor sustituible;
- como primer adaptador futuro, una sola integración HTTP no streaming con el formato OpenAI Chat Completions y respuesta estructurada JSON;
- un contrato de solicitud y respuesta mínimo cuyo único campo molecular generado requerido es `smiles`;
- validación mediante la ruta SMILES aislada ya existente en ChemIO y resultados estables;
- límites de seguridad, privacidad, observabilidad, arquitectura y pruebas offline con proveedores falsos.

En esta sesión sólo se crean y validan artefactos OpenSpec. No se añade código de producto, pruebas de implementación, módulo al catálogo ni dependencia.

## Fuera de alcance

No se implementan cliente HTTP, llamada real a modelo, configuración/UI de API keys, UI/command de usuario, llama.cpp/llama-server, integración específica de LM Studio, varios proveedores, agentes, tool calling, ejecución de código, modificación directa del canvas o de `MolGraph`, memoria, conversación persistente, RAG, búsqueda web, nomenclatura química completa, generación de reacciones, edición molecular iterativa, feedback/optimización Clean2D, multimodalidad, integración con ChemName, benchmarking, cambios grandes de UI ni Campaign 5–9 de Clean2D. Tampoco se modifican algoritmos, scoring, candidatos, políticas o corpus/regresiones de Clean2D.

## Módulos probables

- **Nuevo módulo propuesto, M23 `molecular_assistant` (`src/chemuson/molecular_assistant/`):** dueño de la solicitud de aplicación, abstracción/adaptador de proveedor, decodificación estricta y resultado trazable. Tanto el nombre como la ruta son una propuesta para Phase 1, no un paquete creado en esta sesión. No importa GUI ni Clean2D.
- **M01 `chemio`:** mantiene el parsing/validación de SMILES existente; no adquiere dependencia de IA.
- **M00 `core`:** tipo de salida `MolGraph`.
- **M08/M10/M09 (`gui`, controllers, canvas):** consumidores futuros, no parte de Phase 1.
- **M02 `clean2d`:** consumidor opcional downstream sólo por una orquestación futura; nunca proveedor ni dependencia de IA.

El catálogo `architecture/modules.yml` no se cambia en este OpenSpec. La futura implementación que establezca M23 deberá actualizar el catálogo, APIs/dependencias y contratos de límites en el mismo cambio que introduzca el paquete.

## Suposición arquitectónica comprobada

El repositorio ya expone `chemuson.chemio.rdkit_io.smiles_to_molgraph`, que intenta el worker RDKit aislado y tiene un fallback existente a RDKit directo. La UI de importación SMILES actual no es sólo ese parser: `TemplateController._import_smiles_graph` intenta primero selección de depiction de M02 (`smiles_to_molgraph_best_depiction_with_report`) y después usa el importador ChemIO; la inserción visual llama `canvas._insert_molgraph`. El servicio IA futuro validará con ChemIO, no llamará la ruta de depiction Clean2D. La integración gráfica posterior deberá decidir explícitamente cómo reutilizar la inserción normal y si solicita Clean2D de manera opcional.

También existe M16 `name2structure`, un resolver de nombres comunes/PubChem con conectores propios, no un asistente lingüístico ni una interfaz de proveedor LLM. No se reutiliza ni se amplía en esta campaña.

## Criterio de éxito documental

El contrato deja definidos los límites y resultados de Phase 1, el uso del parser existente, los controles para salida no confiable, la observabilidad sin persistencia sensible, las dependencias prohibidas y una matriz de pruebas deterministas. La validación estricta OpenSpec pasa sin tocar el catálogo ni el código de ChemUSON.