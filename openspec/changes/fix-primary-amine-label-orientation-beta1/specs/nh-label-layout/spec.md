# Especificación delta — layout de etiquetas N–H (aminas ordinarias)

Capability afectada: `nh-label-layout` (creada por
`fix-amide-urea-nh-label-layout-beta1`, aún no archivada). Este delta refina la
elección de posición del H implícito según el contexto químico del nitrógeno.

## ADDED Requirements

### Requirement: Contexto del nitrógeno gobierna el apilado vertical
El canvas SHALL restringir los candidatos de posición del glifo H implícito de un N
neutro dividido (`NH`, `HN`, `NH2`, `H2N`) según el contexto devuelto por
`_nitrogen_geometry_context`, que depende de conectividad y órdenes de enlace, nunca
de proximidad geométrica: `amine` SHALL limitarse a posiciones horizontales
(`left`/`right`); `amide`, `urea`, `aromatic` y los demás contextos conservan las
cuatro direcciones de `fix-amide-urea-nh-label-layout-beta1`. Amonio/N cargado, N
aromático de anillo y centros de coordinación mantienen su tratamiento previo.

#### Scenario: Reproducción 1,1-diamino-3-fenilpropan-2-ona
- **GIVEN** dos grupos amino unidos al mismo C alfa adyacente a `C(=O)` con enlaces
  inclinados ±75°
- **WHEN** se refrescan las etiquetas de ambos N
- **THEN** ambos muestran `NH₂`/`H₂N` horizontales con subíndice 2
- **AND** ningún H queda apilado verticalmente sobre o bajo su N
- **AND** las longitudes de todos los enlaces siguen siendo 40 unidades centro a centro

#### Scenario: Amina primaria alifática, anilina y amina secundaria
- **GIVEN** etilamina, anilina (N simple sobre carbono aromático, enlace no aromático)
  o dietilamina dibujadas en el canvas
- **WHEN** se calcula la posición del H implícito
- **THEN** el contexto del N es `amine`
- **AND** la posición elegida es `left` o `right`, con el N cerca del extremo del
  enlace y el glifo H en la misma línea horizontal que el N

#### Scenario: Controles amida y urea conservan el contrato
- **GIVEN** acetamida o N,N′-dietilurea con N unido directamente al carbonilo
- **WHEN** se calcula la posición del H implícito
- **THEN** el contexto es `amide`/`urea`
- **AND** las posiciones vertical/esquinas elegidas coinciden con las de
  `fix-amide-urea-nh-label-layout-beta1`

### Requirement: Refresco del contexto vecino al mutar enlaces
Al crear, actualizar el orden o eliminar un enlace, el canvas SHALL refrescar la
reposición del H implícito de los N adyacentes a los extremos de ese enlace, además
de los propios extremos, para que la clasificación por contexto nunca quede obsoleta.

#### Scenario: Carbonilo dibujado después del enlace N–C
- **GIVEN** un N ya unido a un carbono al que luego se le añade `C=O`
- **WHEN** se crea el enlace doble al oxígeno
- **THEN** la etiqueta del N recalcula su H con el contexto actualizado
- **AND** el trazado del enlace no atraviesa N ni H y conserva el despeje previo

#### Scenario: Edición, arrastre, undo/redo y tema
- **GIVEN** una amina con H implícito horizontal
- **WHEN** se añade/elimina un enlace, se arrastra el átomo, se deshace/rehace o se
  cambia el tema claro/oscuro
- **THEN** la posición recalculada sigue siendo horizontal y sin solapes
- **AND** el grafo, las cargas y los conteos de H no cambian
