# M01 - Chemio

## Responsabilidad
ChemIO gestiona formatos químicos, parsing, exportación, persistencia y workers RDKit aislados. No contiene scoring, ranking ni layout 2D.

## Dependencias
- Runtime: `core` (M00).
- TYPE_CHECKING: ninguna dependencia hacia otros módulos ChemUSON.
- No depende de Clean2D.
- No depende de GUI ni de Canvas (M09).

## Persistencia GUI-neutral

`persistence.py` define el Protocol interno `PersistenceDocument`, limitado a
`model`, `get_persistence_data`, `load_persistence_data` y
`rebuild_persistence_view`. `PersistenceManager` opera exclusivamente sobre ese
contrato estructural. En runtime puede recibir un `ChemusonCanvas`, pero ChemIO
no importa, nombra ni hereda una clase GUI concreta.

El catálogo de M01 contiene sólo M00 en `current_dependencies` y
`target_dependencies`, sin excepciones temporales ni ciclos.
