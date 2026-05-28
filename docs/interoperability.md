# Interoperabilidad de Formatos

Esta matriz describe el soporte actual de intercambio químico. Los formatos marcados como parciales degradan de forma segura y no deben asumirse como round-trip industrial completo.

| Formato | Importar | Exportar/copiar | Estado | Notas |
| --- | --- | --- | --- | --- |
| `.cmsn` | Sí | Sí | Nativo | Conserva modelo, anotaciones, estilos y semántica Chemuson. |
| SMILES | Sí | Sí | Parcial | Usa RDKit aislado cuando está disponible y fallback interno para subconjuntos. No conserva layout ni anotaciones. |
| Molfile `.mol` | Sí | Sí | Base | Conserva conectividad, coordenadas 2D y parte de atributos químicos. Algunas extensiones visuales degradan. |
| SDF `.sdf` | Sí | Sí | Base/parcial | Se trata como Molfile para estructura principal; propiedades SDF avanzadas no se modelan todavía. |
| CML `.cml` | Sí | Sí | Inicial semántico | Disponible en abrir, guardar, exportar y copiar. Conserva átomos, enlaces, coordenadas 2D, carga, isótopos, radicales y extensiones Chemuson básicas. |
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
