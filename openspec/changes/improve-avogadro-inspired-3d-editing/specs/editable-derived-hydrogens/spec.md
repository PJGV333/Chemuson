## ADDED Requirements

### Requirement: La identidad derivada de H no se confunde con la identidad documental

Los H generados por RDKit SHALL conservar un identificador visual derivado y una asociación inequívoca con su padre. Picking SHALL poder seleccionar y resaltar un H derivado sin tratar su string visual como ID de `MolGraph`. La vista previa de arrastre SHALL modificar únicamente la geometría efímera. Cancelación SHALL dejar intactos el grafo y el historial.

#### Scenario: Selección y arrastre cancelable
- **GIVEN** una escena con un H generado por RDKit
- **WHEN** el usuario lo selecciona y arrastra
- **THEN** el H queda resaltado y su preview cambia en la escena
- **AND** Escape restaura la posición original sin crear un átomo documental ni una entrada undo

### Requirement: Confirmar un H derivado lo promueve mediante una transacción química

Al confirmar una edición persistente de H derivado, ChemUSON SHALL crear un nuevo átomo documental H, enlace con su padre y coordenada 3D en una única operación lógica undoable. SHALL retirar de la metadata efímera ese H y enlace, conservar los otros H derivados y no duplicar el H en generación/optimización RDKit. Undo SHALL restaurar el estado derivado y redo SHALL restaurar la promoción con IDs estables.

#### Scenario: Promoción y ciclo undo/redo
- **GIVEN** un H derivado seleccionado y arrastrado a una posición válida
- **WHEN** el usuario confirma la edición
- **THEN** `MolGraph` contiene un H explícito enlazado al padre con la coordenada confirmada
- **AND** la fórmula/conectividad y escena no contienen un H duplicado
- **WHEN** se hace undo y redo
- **THEN** vuelve primero el H derivado y luego exactamente un H explícito

#### Scenario: Persistencia y Live MM
- **GIVEN** un H promovido con geometría confirmada
- **WHEN** se guarda/reabre `.cmsn` y luego se inicia Live MM
- **THEN** la identidad documental persiste y RDKit conserva su coordenada inicial sin añadir duplicado

### Requirement: Borrar un H derivado no simula un cambio químico

Hasta que una operación validada de deshidrogenación/valencia esté disponible, una solicitud de borrar un H derivado SHALL ser rechazada con un estado claro, sin ocultarlo como si la química hubiese cambiado, sin mutar `MolGraph` y sin añadir undo. H explícitos SHALL seguir el flujo de eliminación química documental existente.

#### Scenario: Borrado de H derivado no soportado
- **GIVEN** un H presente solo en metadata de conformero
- **WHEN** el usuario solicita eliminarlo
- **THEN** la operación se rechaza con explicación de que es derivado y no se altera la escena/grafo de forma engañosa