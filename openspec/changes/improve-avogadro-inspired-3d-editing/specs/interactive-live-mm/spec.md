## ADDED Requirements

### Requirement: Live MM acompaña gestos sin bloquear ni confirmar previews

Con Live MM activo durante construcción/manipulación, ChemUSON SHALL ejecutar RDKit fuera del hilo GUI usando el worker aislado y streaming existente. Los átomos bajo control SHALL mantenerse temporalmente restringidos; los previews SHALL conservar coordenadas/identidades de H y estereoquímica soportada. Nuevas operaciones SHALL coalescer/cancelar solicitudes previas y frames SHALL aplicarse solo si documento, firma química y generación siguen vigentes. Al terminar, una acción lógica SHALL producir como máximo una transacción de coordenadas undoable; frames intermedios no SHALL llenar el historial.

#### Scenario: Construcción y relajación en paralelo
- **GIVEN** una molécula en construcción con Live MM activo
- **WHEN** se añade un átomo/enlace y el usuario continúa manipulando
- **THEN** se observan previews progresivos mientras la UI conserva capacidad de respuesta
- **AND** el fragmento manipulado queda fijo hasta acabar el gesto

#### Scenario: Un trabajo reemplaza a otro
- **GIVEN** un job de optimización activa
- **WHEN** se inicia otra operación o cambia la identidad química
- **THEN** se cancela/desautoriza el job anterior y ninguno de sus frames puede modificar la molécula nueva

#### Scenario: Un commit por gesto
- **GIVEN** una secuencia de frames de Live MM para un único gesto
- **WHEN** termina y se confirma la relajación vigente
- **THEN** undo restaura la geometría previa completa y redo la geometría final
- **AND** frames intermedios no aparecen como entradas independientes