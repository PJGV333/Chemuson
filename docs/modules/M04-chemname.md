# M04: chemname

**Responsabilidad**: Motor de nomenclatura IUPAC-lite. Resuelve nombres químicos para un subconjunto de estructuras (hidrocarburos, aromáticos sencillos, heterociclos comunes) a partir de la topología del grafo.

**APIs Principales**:
- `chemuson.chemname.engine`: Función principal `iupac_name`.
- `chemuson.chemname.options`: Configuración de la nomenclatura (`NameOptions`).

**Riesgos**:
- **Ciclo con `chemcalc`**: Existe una circularidad con `chemcalc` (debido a que la nomenclatura requiere el recuento de hidrógenos implícitos de `valence.py`).
- **Alcance limitado**: Actualmente solo cubre un subconjunto de la química orgánica (IUPAC-lite); no es un motor de nomenclatura completo.
