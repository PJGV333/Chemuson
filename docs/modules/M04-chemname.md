# M04 - ChemName

## Responsabilidad
`chemname` implementa el sistema de nomenclatura IUPAC-lite de Chemuson. Proporciona una solución semisistemática para la identificación de estructuras, incluyendo la selección de cadenas principales, la identificación de anillos, el manejo de grupos funcionales, la determinación de la estereoquímica y la gestión de la coordinación.

## APIs Principales
- **Nomenclatura**: `iupac_name`.
- **Configuración**: `NameOptions`.

## Riesgos Conocidos
La cobertura de este módulo es parcial y está diseñada para ser un soporte ágil. Cualquier modificación en la lógica de selección de cadenas o en la identificación de grupos funcionales debe ser validada con extremo cuidado para no romper la cobertura histórica de los tests de regresión.

## Dependencias
- Runtime: `core` (M00), `chemcalc` (M03), `chemio` (M01, lazy), `utils` (M15).
- Ciclo conocido: `chemcalc` ↔ `chemname` vía `chemcalc` que importa `MolView` de `chemname`, y `chemname` que importa `implicit_h_count` de `chemcalc`.
