# Spec Delta — interoperabilidad de archivos MOL

## ADDED Requirements

### Requirement: Guardar como MOL conserva el filtro y el sufijo
Cuando se elige el filtro MOL, el sistema SHALL escribir un archivo `.mol`, añadiendo el sufijo si no existe y evitando duplicarlo. La ruta del documento, estado limpio, recientes y limpieza de autosave SHALL actualizarse sólo tras escritura correcta.

#### Scenario: Guardado MOL sin sufijo
- **WHEN** el usuario selecciona el filtro MOL y escribe `molecula`
- **THEN** el archivo se guarda como `molecula.mol` y el documento conserva esa ruta

#### Scenario: Escritura fallida
- **WHEN** el escritor o el sistema de archivos rechaza el guardado MOL
- **THEN** se muestra el error concreto y el documento sigue sin guardar, con su ruta previa

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

- `.cmsn` y sus payloads no cambian.
- Ninguna exportación cambia el grafo/documento fuente.
- Sin reflexión ni asignación estéreo ausente.
- No cambian Clean2D, ChemName, reglas IUPAC ni dependencias runtime.
