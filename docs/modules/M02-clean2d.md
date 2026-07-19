# M02 - Clean2D

## Responsabilidad
`clean2d` proporciona el motor de limpieza y representación (depiction) 2D de Chemuson. El módulo está diseñado para ser desacoplado de la GUI, operando sobre modelos de grafo y coordenadas. Sus funciones incluyen la generación de candidatos de limpieza, el ranking de calidad, la aplicación de invariantes geométricos, la reparación de errores locales y la gestión de políticas para moléculas con alta complejidad estructural.

## APIs Principales
- **Motor de Limpieza**: `run_clean2d_engine`.
- **Generación de Candidatos**: `generate_clean2d_candidates`.
- **Evaluación y Ranking**: `rank_clean2d_candidates`, `evaluate_clean2d_layout`.
- **Seguridad**: `is_clean2d_candidate_safe`.

## Riesgos Conocidos
Es uno de los módulos más complejos algorítmicamente. Los cambios en las heurísticas de limpieza o en las reglas de ranking pueden alterar la consistencia visual de las moléculas y provocar regresiones en la calidad de la representación, lo cual es crítico para la utilidad científica de la herramienta.

## Dependencias
- Runtime: `core` (M00), `chemio` (M01) — las dependencias de chemio son imports lazy dentro de funciones (`rdkit_safe`), no en nivel de módulo.
- Ciclo conocido: `chemio` ↔ `clean2d` vía `depiction_candidates.py` (chemio) y `block_unwrap.py`/`scaffold_depiction.py` (clean2d).
