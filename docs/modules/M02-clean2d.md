# M02: clean2d

**Responsabilidad**: Servicios de limpieza y optimización de disposiciones 2D (depiction) desacoplados de la GUI. Incluye motores de búsqueda de candidatos de layout, validación de calidad (ciclos, distancias, degeneración) y reporte de defectos locales.

**APIs Principales**:
- `chemuson.clean2d.v2`: Parámetros y optimización.
- `chemuson.clean2d.safety`: Validación de calidad y seguridad geométrica.
- `chemuson.clean2d.engine`: Motor de generación de candidatos y ejecución de procesos de limpieza.
- `chemuson.clean2d.quality_reporting`: Reportes de diagnóstico y métricas de calidad.

**Riesgos**:
- **Dependencia de RDKit**: Utiliza RDKit para ciertas operaciones de limpieza y validación (ej. `run_clean2d_engine`).
- **Ciclo con `chemio`**: Existe una circularidad con `chemio` (debido a `depiction_candidates.py`) que debe ser resuelta.
- **Complejidad algorítmica**: Los motores de búsqueda de candidatos y optimización pueden ser costosos computacionalmente.
