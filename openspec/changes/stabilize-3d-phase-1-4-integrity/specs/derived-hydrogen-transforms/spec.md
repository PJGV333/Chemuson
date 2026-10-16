## ADDED Requirements

### Requirement: Los hidrógenos derivados acompañan transformaciones rígidas

Las coordenadas de H generados por RDKit SHALL permanecer en metadata efímera enlazada a su padre documental. Cuando una operación aplique un traslado o una rotación rígida a un fragmento documental, SHALL aplicar exactamente la misma transformación a cada H derivado enlazado a un átomo de ese fragmento. SHALL conservar distancias internas C–H, direcciones de enlaces y ángulos del fragmento bajo transformaciones rígidas. Una edición local no rígida SHALL conservar el H unido a su padre y SHALL indicar/evitar afirmar una reorientación no calculada.

#### Scenario: Rotar fragmento completo
- **GIVEN** un conformero de etano, etanol o butano con H derivados
- **AND** un fragmento rígido que incluye uno o más padres de H
- **WHEN** el fragmento rota alrededor de un eje 3D
- **THEN** sus H reciben la misma transformación de Rodrigues que el fragmento
- **AND** las distancias padre–H y los ángulos internos permanecen invariantes
- **AND** los vectores C–H rotan respecto al sistema global

#### Scenario: Torsionar butano
- **GIVEN** butano con H derivados en ambos lados de un enlace central acíclico
- **WHEN** se cambia su dihedro
- **THEN** los H del lado móvil rotan con ese lado y los del lado fijo no se mueven
- **AND** no se altera `MolGraph` ni se crean IDs documentales para esos H

#### Scenario: Undo, cancelación y relajación posterior
- **GIVEN** una transformación que mueve padres e H derivados
- **WHEN** se deshace, rehace, cancela o inicia una relajación Live MM válida
- **THEN** la escena reconstruida mantiene la asociación padre–H y el estado comprometido es consistente
- **AND** el archivo `.cmsn` contiene sólo átomos documentales

#### Scenario: Cambio de conectividad
- **GIVEN** metadata de H derivados cuya adyacencia ya no coincide con la topología
- **WHEN** cambia la conectividad documental
- **THEN** los H incompatibles no se reutilizan y pueden regenerarse desde RDKit