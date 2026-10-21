# Propuesta — tipografía de etiquetas de diagramas de energía

## Problema

Las opciones visibles de presets electrónicos presentan fórmulas diatómicas y configuraciones electrónicas con números/cargas en línea base; algunas etiquetas españolas omiten tildes. Esto dificulta distinguir, por ejemplo, superóxido (`O₂⁻`) de peróxido (`O₂²⁻`).

## Alcance

- Corregir exclusivamente el texto visible de presets en la paleta y submenús de diagramas de energía, incluyendo tooltips pertinentes.
- Presentar `B₂`, `C₂`, `N₂`, `O₂`, `F₂`, `O₂⁻`, `O₂²⁻`, `Híbrido sp`, `Híbrido sp²`, `Híbrido sp³` y las opciones `d⁴`, `d⁵`, `d⁶`, `d⁸` conservando geometría y espín.
- Mantener `2s`, `2p` y otros números cuánticos principales en línea base.
- Añadir pruebas de identidad de preset, presentación, equivalencia científica, selección Qt y persistencia; capturas comparativas antes/después en tema claro/oscuro y escala disponible.

## Invariantes y fuera de alcance

No cambiar IDs como `mo:O2^2-`, `hybrid_sp3` o `ligand_field:d6_octahedral_low`; valores ni propiedades científicas de presets; electrones, cargas, ocupaciones, niveles, energías, degeneración, orden orbital o geometría; construcción/edición de diagramas; ni el formato/títulos persistidos de documentos. No reformar fuentes, colores, dimensiones o el tema general. Sin dependencias nuevas, cambios de paquetes ni funcionalidades nuevas.
