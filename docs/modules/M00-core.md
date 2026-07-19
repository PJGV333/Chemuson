# M00 - Core

## Responsabilidad
El módulo `core` es el corazón de Chemuson. Proporciona el modelo molecular fundamental, el estado químico (`ChemState`), las estructuras de datos para el grafo químico multicapa (`MolGraph`, `BlockGraph`, etc.) y las herramientas esenciales de análisis elemental. Es un subsistema puramente lógico, sin dependencias de la interfaz de usuario (GUI) ni de motores de dibujo externos.

## APIs Principales
- **Modelado**: `Atom`, `Bond`, `BondStyle`, `BondStereo`, `ChemState`, `MolGraph`.
- **Grafos Multicapa**: `BlockGraph`, `InteractionGraph`, `MotifGraph`, `build_multilayer_chemical_graph`.
- **Análisis Elemental**: `elemental_percentages`, `parse_formula`, `molecular_weight`, `find_solvate_candidates`.

## Riesgos Conocidos
Al ser la base sobre la que se construyen todos los demás paquetes, cualquier cambio en la definición de un átomo, enlace o en la estructura del grafo tiene un impacto sistémico. Se requiere una cobertura de tests extremadamente rigurosa para evitar regresiones en la integridad del modelo químico.

## Dependencias
Core no importa ningún otro módulo chemuson en runtime. Es el único módulo con `current_dependencies: []`.
