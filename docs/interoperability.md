# Interoperabilidad de Formatos

Esta matriz describe el soporte actual de intercambio químico. Los formatos marcados como parciales degradan de forma segura y no deben asumirse como round-trip industrial completo.

| Formato | Importar | Exportar/copiar | Estado | Notas |
| --- | --- | --- | --- | --- |
| `.cmsn` | Sí | Sí | Nativo | Conserva modelo, anotaciones, estilos y semántica Chemuson. |
| SMILES | Sí | Sí | Parcial | Usa RDKit aislado cuando está disponible y fallback interno para subconjuntos. No conserva layout ni anotaciones. |
| Molfile `.mol` | Sí | Sí | Base | Conserva conectividad, coordenadas 2D y parte de atributos químicos. Algunas extensiones visuales degradan. |
| SDF `.sdf` | Sí | Sí | Base/parcial | Se trata como Molfile para estructura principal; propiedades SDF avanzadas no se modelan todavía. |
| CML `.cml` | Sí | Sí | Inicial semántico | Disponible en abrir, guardar, exportar y copiar. Conserva átomos, enlaces, coordenadas 2D, carga, isótopos, radicales y extensiones Chemuson básicas. |
| XYZ `.xyz` | No | Sí | 3D/compchem | Exporta el `CoordinateSet3D` activo del dock 3D. No conserva conectividad ni semántica 2D. |
| ORCA/Gaussian/NWChem input | No | Sí | CompChem MVP | Exporta geometría 3D activa, carga y parámetros básicos. No intenta validar el método/base para cada paquete. |
| SVG | No químico | Sí | Publicación | Exporta visual; no es round-trip químico. |
| PDF | No químico | Sí | Publicación | Exporta visual; no es round-trip químico. |
| PNG | No químico | Sí | Publicación | Raster visual; no conserva semántica química. |
| CDXML | No | No | Pendiente | Siguiente fase de interoperabilidad industrial. |
| MRV | No | No | Pendiente | Siguiente fase de interoperabilidad industrial. |
| PDB | No | No | Pendiente | Requiere modelo 3D/biomolecular más completo. |

## CML Inicial

El soporte CML actual usa CML estándar para `atomArray`/`bondArray` y extensiones namespaced `chemuson:*` para datos que no tienen equivalente directo en el subconjunto usado.

Se preserva:
- Elemento, coordenadas `x2`/`y2`, carga formal, isótopos y radicales.
- Orden de enlace y aromaticidad básica (`order="A"`).
- Estilo/estereo básico de Chemuson como metadatos namespaced.
- Sustituyentes Markush de átomos `R/Rn` mediante `chemuson:rGroupSubstituents`.

Limitaciones actuales:
- No exporta todavía flechas, textos, corchetes ni objetos gráficos como CML semántico.
- No intenta interpretar toda la especificación CML ni diccionarios externos.
- CDXML/MRV siguen pendientes para compatibilidad con suites comerciales.

## Limitaciones explícitas de junio 2026

- CDXML, MRV y PDB siguen sin soporte.
- Los espectros son heurísticos/estimados y no forman parte de un formato experimental validado.
- IUPAC-lite es parcial; ante estructuras fuera de alcance degrada a `N/D`.
- Las exportaciones compchem dependen de que exista un conformero 3D activo; Open Babel es opcional y requiere `obabel` en `PATH`.
