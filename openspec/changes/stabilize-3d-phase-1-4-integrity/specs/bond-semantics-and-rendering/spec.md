## ADDED Requirements

### Requirement: Aromaticidad y orden de enlace permanecen semánticamente consistentes

La conversión ChemUSON→RDKit SHALL preservar aromaticidad explícita, orden de enlace y conectividad. Una forma Kekulé válida SHALL seguir siendo equivalente para cálculo químico después de sanitización RDKit; SHALL NOT marcarse aromático un ciclo sólo por ser conjugado. La escena 3D SHALL retener orden y bandera aromática para renderizar enlaces simples, dobles, triples y aromáticos de manera distinguible sin mutar el documento.

#### Scenario: Aromático explícito y Kekulé equivalente
- **GIVEN** benceno con enlaces marcados aromáticos y benceno Kekulé alternante
- **WHEN** cada grafo se convierte y sanitiza en RDKit
- **THEN** ambos se reconocen como el mismo sistema aromático
- **AND** se conserva la distinción documental original y el grafo no muta

#### Scenario: Controles no aromáticos
- **GIVEN** ciclohexano, eteno y etino
- **WHEN** se convierte, conforma y construye la escena
- **THEN** ciclohexano queda saturado, eteno conserva un doble enlace y etino un triple enlace
- **AND** el render de escena representa la multiplicidad correspondiente

#### Scenario: Enlaces aromáticos explícitos
- **GIVEN** un anillo con bandera aromática explícita y orden Kekulé/documental almacenado
- **WHEN** se construye la escena y la malla
- **THEN** la aromaticidad no se pierde en el modelo de escena y el anillo se representa con una convención uniforme
- **AND** picking, conectividad y persistencia documental no cambian