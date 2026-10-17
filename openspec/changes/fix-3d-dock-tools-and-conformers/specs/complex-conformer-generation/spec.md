## ADDED Requirements

### Requirement: Generación 3D de moléculas complejas es acotada y reproducible

La generación SHALL continuar en worker aislado fuera del hilo GUI, con timeout configurable y cancelación efectiva del proceso. Los intentos ETKDG SHALL ser finitos, diferenciados y trazables; opciones de macro ciclo SHALL usarse solo cuando corresponda. Geometría aceptada SHALL tener coordenadas finitas y cumplir estereoquímica documentada. Si el embedding válido sobrevive pero MMFF/UFF no tiene parámetros o falla, SHALL poder publicarse como `unoptimized` con explicación clara.

#### Scenario: Embedding válido y minimización indisponible
- **GIVEN** un grafo que permite embedding pero no cuenta con parámetros MMFF/UFF
- **WHEN** termina el intento limitado de generación
- **THEN** se presenta una geometría 3D válida marcada como no optimizada y el motivo del fallo de minimización

#### Scenario: Timeout/cancelación de estructura difícil
- **GIVEN** un cálculo de conformero asíncrono activo
- **WHEN** vence el timeout o el usuario cancela
- **THEN** el proceso RDKit termina, se informa `timeout`/`cancelled` con etapa alcanzada y la GUI sigue disponible

### Requirement: Fallos no se ocultan como viewport vacío

Cada resultado SHALL informar etapa, conteos del grafo, cargas/valencias, estereoquímica, tiempos disponibles, retorno RDKit, timeout/cancelación, estado snapshot y estado de escena/render. Los fallos SHALL diferenciar al menos `invalid_graph`, `embed_failed`, `timeout`, `unsupported_stereochemistry`, ausencia de parámetros, worker/scene/display y resultado obsoleto. Viewport y panel CompChem SHALL mostrar mensaje legible y opción de reintento cuando proceda; resultados obsoletos SHALL nunca sobrescribir el documento activo.

#### Scenario: Fallo de embedding
- **GIVEN** una estructura no embebible
- **WHEN** finalizan los intentos finitos
- **THEN** el panel y viewport muestran `embed_failed`, código/motivo de RDKit, etapa y resumen químico, no solo “Sin geometría 3D”

#### Scenario: Trabajo obsoleto
- **GIVEN** una generación en pestaña A
- **WHEN** cambia la topología o se activa la pestaña B antes de la respuesta
- **THEN** el resultado no se aplica a A/B y el estado diagnóstico identifica que fue descartado

### Requirement: El progreso representa etapas reales y la escena se verifica por capas

El usuario SHALL recibir transiciones reales (payload/sanitización, H/estereoquímica, embedding, minimización, publicación) sin porcentajes ficticios. La validación SHALL separar generación de CoordinateSet, construcción de SceneMolecule3D, intento de upload/render y aceptación visual GPU. El caso autenticado tetrandrina y otros grafos deben cubrir el camino sin depender del archivo local del usuario.

#### Scenario: Tetrandrina de referencia
- **GIVEN** el isomeric SMILES real documentado para PubChem CID 73078
- **WHEN** se ejecuta el worker aislado y la publicación a escena
- **THEN** se verifica el número de átomos/enlaces, estereoquímica, coordenadas finitas y escena
- **AND** se declara por separado que no se validó el `.cmsn` original ni el render real si no se inició GUI
