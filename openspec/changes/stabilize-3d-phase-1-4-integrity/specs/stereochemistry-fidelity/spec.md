## ADDED Requirements

### Requirement: La generación 3D respeta estereoquímica explícita inequívoca

El worker RDKit SHALL recibir la estereoquímica documental disponible para centros tetraédricos y dobles enlaces. La importación SMILES SHALL conservar los descriptores explícitos R/S y E/Z en `MolGraph`, mapeándolos desde índices RDKit a IDs estables; no SHALL degradarlos a una cuña/hash visual ambigua ni descartarlos. Para configuraciones R/S o E/Z explícitas que puedan codificarse sin ambigüedad, el worker SHALL asignar la representación RDKit correspondiente y verificarla después de conformar y optimizar. SHALL NOT aceptar silenciosamente una configuración opuesta. Si cuñas/hash, CIP u otro dato son insuficientes o contradictorios para reconstruir la configuración, el servicio SHALL devolver un diagnóstico explícito y no etiquetar el conformero como estereoquímicamente fiel.

#### Scenario: Centros R y S
- **GIVEN** dos grafos con configuración tetraédrica inequívoca y CIP objetivo R o S
- **WHEN** se construye el payload RDKit y se genera un conformero
- **THEN** el CIP resultante coincide con el objetivo documentado
- **AND** el grafo original conserva sus IDs, coordenadas 2D y etiquetas

#### Scenario: Doble enlace E/Z
- **GIVEN** un doble enlace con dos sustituyentes distinguibles y configuración E o Z explícita
- **WHEN** el worker reconstruye y optimiza la molécula
- **THEN** RDKit conserva la configuración objetivo y ChemUSON la verifica antes y después de la minimización

#### Scenario: Sin especificación o especificación ambigua
- **GIVEN** una molécula sin configuración especificada, con sustituyentes equivalentes o con cuña/hash insuficientes
- **WHEN** se solicita un conformero
- **THEN** la configuración no se inventa ni se declara preservada
- **AND** la ambigüedad se conserva como no especificada o se devuelve un diagnóstico controlado