# User Manual Documentation Specification

## Purpose

Garantizar que el manual de usuario de ChemUSON describa únicamente funciones verificadas de la versión documentada, con instrucciones reproducibles, atajos correctos e ilustraciones auténticas.

## ADDED Requirements

### Requirement: Verificación funcional de cada afirmación
El manual SHALL documentar solo funciones, rutas de menú, atajos y comportamientos comprobados contra el código y la ejecución de la versión de referencia, y SHALL indicar explícitamente cuando una función sea experimental, esté incompleta o dependa de componentes externos.

#### Scenario: Función no verificable
- **GIVEN** una función mencionada en una revisión anterior del manual
- **WHEN** no puede comprobarse en la versión de referencia
- **THEN** el manual no la describe como disponible y el hueco queda registrado en la evidencia del cambio

#### Scenario: Función experimental
- **GIVEN** una función presente en la interfaz pero que requiere configuración externa
- **WHEN** el usuario la busca en el manual
- **THEN** el manual la describe marcándola como experimental y explica qué configuración necesita

### Requirement: Un único desarrollo por función
El manual SHALL explicar cada función en un solo lugar, con su finalidad, su ubicación, su procedimiento, un ejemplo y sus limitaciones, y SHALL enlazar el resto de referencias mediante anclas internas.

#### Scenario: Función citada en varios contextos
- **GIVEN** una función usada desde menús, rail y menú contextual
- **WHEN** el manual la menciona en más de un capítulo
- **THEN** solo un capítulo desarrolla el procedimiento y los demás enlazan a él

### Requirement: Atajos coherentes con la interfaz
El manual SHALL listar únicamente atajos que la interfaz de la versión de referencia muestre o procese, y SHALL advertir cuando una entrada de menú no tenga atajo visible.

#### Scenario: Atajo declarado en versiones anteriores del manual
- **GIVEN** un atajo documentado en una revisión anterior
- **WHEN** la versión de referencia no lo asocia a esa acción
- **THEN** el manual lo retira y lo señala como no asignado

### Requirement: Ilustraciones auténticas y trazables
Toda figura del manual SHALL ser una captura de la aplicación de la versión de referencia o una ilustración técnica del proyecto identificada como tal, almacenada en una carpeta ordenada del repositorio y documentada con su procedencia y su método de regeneración.

#### Scenario: Captura no disponible
- **GIVEN** una función cuya interfaz no puede capturarse de forma fiable
- **WHEN** se redacta su capítulo
- **THEN** el capítulo se publica sin figura y declara que no existe captura auténtica

#### Scenario: Evidencia interna de desarrollo
- **GIVEN** una galería generada por herramientas internas de comparación
- **WHEN** no aporta valor didáctico al usuario
- **THEN** deja de enlazarse desde el manual sin borrar la evidencia del repositorio

### Requirement: Separación de manipulación 2D, fragmentos y modelo 3D
El manual SHALL distinguir entre giro del dibujo en el plano, reorientación de ramas o fragmentos, proyección de perspectiva tridimensional sobre un esquema plano y modelo tridimensional real.

#### Scenario: Usuario busca girar una molécula
- **GIVEN** una estructura dibujada en el lienzo
- **WHEN** consulta el capítulo de rotación
- **THEN** encuentra una tabla que asocia cada objetivo con la entrada correcta y su capítulo
