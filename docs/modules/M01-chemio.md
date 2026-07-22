# M01 - Chemio

## Responsabilidad
ChemIO gestiona formatos químicos, parsing, exportación, persistencia y workers RDKit aislados. No contiene scoring, ranking ni layout 2D.

## Dependencias
- Runtime: `core` (M00).
- TYPE_CHECKING: `gui.canvas` (M09), exclusivamente en `persistence.py`.
- No depende de Clean2D.
