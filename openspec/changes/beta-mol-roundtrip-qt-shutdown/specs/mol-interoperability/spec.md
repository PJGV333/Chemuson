# Spec Delta — interoperabilidad de archivos MOL

## ADDED Requirements

### Requirement: Exportar MOL no guarda el documento completo
La acción de Guardar SHALL persistir el documento completo en `.cmsn`. La exportación `.mol` SHALL presentarse por separado como exportación únicamente de la estructura molecular. El archivo SHALL contener sólo un registro CTAB; texto y flechas no están representados y la interfaz MUST indicarlo explícitamente. Exportar SHALL añadir `.mol` idempotentemente y MUST NOT cambiar la ruta nativa, clean/dirty del undo stack, estado/origen/temporizadores de autosave ni recientes. Un intento de tratar `.mol` como destino de Guardar SHALL ser rechazado con una explicación y no SHALL alterar el estado del documento.

#### Scenario: Documento con molécula, texto y flecha
- **WHEN** el usuario exporta ese documento como MOL
- **THEN** el CTAB contiene sólo la molécula; se informa que texto/flechas no se incluyen y que `.cmsn` conserva el documento completo
- **AND** ruta nativa, dirty state, autosave y recientes permanecen intactos

#### Scenario: Exportación MOL sin sufijo o con sufijo
- **WHEN** el usuario exporta como `molecula` o `molecula.mol`
- **THEN** se escribe exactamente `molecula.mol` sin alterar el estado del documento

#### Scenario: Escritura MOL fallida
- **WHEN** el escritor o el sistema de archivos rechaza la exportación
- **THEN** se muestra el error concreto y no cambian ruta, estado dirty, autosave ni recientes

#### Scenario: Usuario intenta Guardar en una ruta MOL
- **WHEN** el destino de Guardar es `.mol`
- **THEN** la acción no escribe el archivo ni marca guardado el documento, e informa de usar Exportar estructura molecular

### Requirement: Coordenadas MOL 2D usan escala química
El exportador SHALL preservar proporciones, ángulos, conectividad y orientación, usando un escalado positivo uniforme que sitúe la mediana de longitud de enlace estructural cerca de 1.5 unidades MOL. La conversión SHALL aplicar la misma regla a todos los backends soportados y no SHALL mutar las coordenadas del documento.

#### Scenario: Exportar coordenadas del lienzo
- **WHEN** un MolGraph en escala de lienzo se exporta como MOL
- **THEN** las longitudes del CTAB tienen escala química convencional y el grafo fuente conserva exactamente sus coordenadas

### Requirement: Round-trip mantiene identidad química y stereo
Para estructuras soportadas, `MolGraph → MOL → MolGraph` SHALL preservar elementos, número de átomos/enlaces, conectividad, órdenes, cargas, isótopos y estereoquímica explícita. Configuración CIP definida SHALL ser equivalente; el enantiómero opuesto SHALL seguir siendo distinto. Centro no especificado SHALL permanecer no asignado. El backend SHALL fallar con diagnóstico explícito si no puede cumplir el contrato.

#### Scenario: Fixture S y control R
- **WHEN** se importan y exportan MOL para `C[C@H](F)CO` y `C[C@@H](F)CO`
- **THEN** se preservan fórmula, grafo molecular y sus configuraciones CIP opuestas

#### Scenario: Stereo ausente
- **WHEN** el MOL contiene un centro potencial sin asignación
- **THEN** ChemIO no inventa configuración ni cuña/hash definido

### Requirement: Apertura desde archivo usa el flujo real
El controlador SHALL abrir archivos MOL generados por ChemUSON y RDKit. El error de lectura SHALL ser visible e identificar la excepción/etapa; no SHALL sustituir ni marcar guardado un documento anterior al fallar.

#### Scenario: Archivos MOL externos e internos
- **WHEN** se abre un archivo MOL válido de RDKit o ChemUSON mediante el controlador
- **THEN** se reconstruye el grafo equivalente en el canvas

### Requirement: SDF no se anuncia parcialmente
Hasta implementar el contrato SDF multi-registro, la GUI SHALL NOT anunciar ni escribir MOL simple como SDF. La ruta `.sdf` explícita SHALL rechazarse con un mensaje claro.

#### Scenario: SDF no implementado
- **WHEN** se ofrece un archivo `.sdf` mediante el flujo ChemUSON
- **THEN** se informa que SDF no está soportado y no se guarda un MOL mal etiquetado

## Invariantes

- `.cmsn` y sus payloads no cambian; Guardar conserva el documento completo.
- MOL se presenta únicamente como exportación molecular parcial, no como copia completa del documento.
- Ninguna exportación cambia el grafo/documento fuente, sus rutas, estado dirty o autosave.
- Sin reflexión ni asignación estéreo ausente.
- No cambian Clean2D, ChemName, reglas IUPAC ni dependencias runtime.
