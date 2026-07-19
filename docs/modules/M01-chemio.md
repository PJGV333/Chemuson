# M01: chemio

**Responsabilidad**: Importación/exportación de datos químicos y persistencia del estado del editor. Maneja la traducción entre el modelo interno de ChemUSON y formatos estándar (SMILES, MOL, SVG) mediante RDKit u otros motores.

**APIs Principales**:
- `chemuson.chemio.rdkit_io`: Conversión de grafos a formatos externos.
- `chemuson.chemio.persistence.PersistenceManager`: Serialización y carga del estado `.cmsn`.
- `chemuson.chemio.depiction_candidates`: Generación de candidatos para visualización.

**Riesgos**:
- **Dependencia de RDKit**: Dependencia externa crítica que se maneja de forma perezosa (lazy) para evitar problemas de compatibilidad en el proceso principal.
- **Ciclo con `clean2d`**: Existe una circularidad conocida con `clean2d` (debido a `depiction_candidates.py`) que debe ser resuelta.
- **Complejidad de persistencia**: La restauración completa del lienzo requiere una sincronización perfecta entre el modelo y la vista (mixins de la GUI).
