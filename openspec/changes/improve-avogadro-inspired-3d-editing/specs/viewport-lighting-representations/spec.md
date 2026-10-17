## ADDED Requirements

### Requirement: La iluminación es estable respecto a la cámara y a la orientación molecular

El renderer SHALL transformar las normales al mismo espacio de coordenadas que luces y cámara, y usar iluminación key/fill/ambient equilibrada con especular moderada. SHALL conservar colores químicos por elemento; no SHALL resolver sombras incorrectas aclarando globalmente los carbonos. La corrección técnica del shader SHALL tener una prueba de contrato y la legibilidad final SHALL identificarse separadamente como validación visual.

#### Scenario: Girar molécula
- **GIVEN** carbono, hidrógeno y enlaces visibles en una escena
- **WHEN** se rota la cámara a orientaciones distintas
- **THEN** la contribución difusa/especular permanece consistente en espacio de vista y los colores elementales siguen reconocibles

### Requirement: Los estilos básicos conservan picking y semántica de enlace

ChemUSON SHALL ofrecer wireframe, bolas y varillas y espacio lleno aproximado con radios documentados. El estilo activo SHALL gobernar tanto malla visible como targets de picking. Órdenes dobles, triples y aromaticidad explícita SHALL conservar su representación semántica cuando el estilo muestre enlaces. Espacio lleno no SHALL presentarse como una superficie electrónica o energética.

#### Scenario: Cambiar estilo sin romper selección
- **GIVEN** una molécula con H, enlaces múltiples y aromaticidad
- **WHEN** el usuario alterna los tres estilos
- **THEN** cada estilo muestra su geometría definida, picking selecciona solo elementos visibles y undo/documento/persistencia química no cambian

#### Scenario: Límite de superficies
- **GIVEN** una vista de espacio lleno
- **WHEN** se inspecciona la escena
- **THEN** la UI la identifica como representación de átomos por radios aproximados y no como densidad electrónica o potencial