# Especificación — layout de etiquetas N–H

## Propósito

Mantener la continuidad visual de enlaces y grupos `N–H` sin modificar el centro del átomo, el modelo químico o la geometría del enlace.

## ADDED Requirements

### Requirement: N anclado y H reposicionable
El canvas SHALL anclar el centro visual del glifo N al centro del átomo cuando su label contiene H inline (`NH`, `HN`, `NH2`, `H2N`). SHALL ubicar el H/H₂ a un lado, arriba o abajo según direcciones de enlaces estructurales y espacio gráfico; las letras no se rotan. La selección no depende de Clean2D ni cambia la etiqueta química, coordenadas o conteos de H.

#### Scenario: Urea N,N'-sustituida
- **GIVEN** la estructura se dibuja con eventos reales de anillo, clics y arrastres del canvas
- **WHEN** se actualizan las etiquetas de ambos N alrededor de `C(=O)`
- **THEN** ambos N permanecen centrados en sus respectivos átomos
- **AND** los H pueden ocupar lados distintos, arriba o abajo según el espacio libre
- **AND** los trazos N–C(=O) no cruzan los glifos y conservan las 40 unidades entre centros

#### Scenario: Amina primaria/secundaria y acetanilida
- **GIVEN** una N de amina o una N de acetanilida tiene H implícito visible en su etiqueta
- **WHEN** vecinos u orientación dejan espacio en otro cuadrante
- **THEN** el H se recoloca independientemente del N sin cambiar valencia ni geometría
- **AND** el label N sin H y las etiquetas O–H/S–H conservan su comportamiento previo

### Requirement: Forma visible íntegra
La forma/bounds del label usada por recorte y colisiones SHALL incluir el H separado; el margen existente al enlace y los tratamientos de estilos no afectados SHALL preservarse.

#### Scenario: Recorte respecto a NH
- **GIVEN** el N y H están renderizados en posiciones separadas
- **WHEN** se calcula el trim del enlace vecino
- **THEN** el path visible no intersecta N, H/H₂, carga ni otro label
- **AND** el trim de la corrección previa no vuelve a la proyección excesiva

### Requirement: Invariantes de dibujo
El layout SHALL preservar coordenadas, distancias/ángulos de enlace, conectividad, estilos, Clean2D, preview/commit, undo/redo y CMSN.

#### Scenario: Edición y persistencia
- **GIVEN** se dibuja un caso con eventos Qt y se captura su grafo
- **WHEN** se confirma, deshace/rehace y guarda/abre
- **THEN** los endpoints de preview/commit y los snapshots del grafo coinciden
- **AND** solo la composición visual de H difiere entre las capturas antes/después
