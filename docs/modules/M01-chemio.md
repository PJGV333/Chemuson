# M01 - Chemio

## Responsabilidad
El módulo `chemio` gestiona la entrada, salida y persistencia de datos químicos. Se encarga de la conversión entre formatos estándar (SMILES, Molfile, CML) y el formato interno de Chemuson (`.cmsn`). Además, actúa como la capa de abstracción segura para interactuar con RDKit mediante trabajadores (workers) que evitan problemas de estabilidad en el proceso principal.

## APIs Principales
- **I/O**: `rdkit_io`, `cml_io`.
- **Persistencia**: `persistence` (gestión de archivos `.cmsn`).
- **Depiction**: `depiction_candidates` (generación de candidatos para limpieza 2D).
- **Seguridad**: `rdkit_safe` (acceso controlado a RDKit).

## Riesgos Conocidos
Depende críticamente de la correcta implementación de los workers para RDKit. Cualquier error en la serialización de la persistencia o en la conversión de formatos puede provocar la pérdida de datos del usuario o la corrupción del estado del editor.

## Dependencias
- Runtime: `core` (M00), `clean2d` (M02).
- TYPE_CHECKING: `gui.canvas` (M09) — solo para tipado en `persistence.py`.
- subprocess: `_rdkit_worker.py` es invocado por `subprocess` desde `rdkit_safe.py`, no aparece en imports normales.
