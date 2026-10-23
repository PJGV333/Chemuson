# Propuesta — orientación horizontal del H implícito en aminas ordinarias

## Problema

En la AppImage integrada de 0.3.0-beta.1, **1,1-diamino-3-fenilpropan-2-ona** muestra
el H₂ de un grupo amino **apilado sobre el N** y el del otro amino **apilado bajo el N**
(es el mismo átomo C1, con enlaces simétricos ±75°). El cambio
`fix-amide-urea-nh-label-layout-beta1` habilitó la estrategia de cuatro direcciones
(izquierda/derecha/arriba/abajo) para **todo** N neutro con etiqueta `NH/HN/NH2/H2N`,
incluidas las aminas ordinarias, y la jerarquía de colisiones elige arriba/abajo cuando
el sector libre lo favorece. El contrato existente en
`tests/test_amide_urea_2d_continuity.py` incluso codificaba `above` para una amina
secundaria. Las aminas ordinarias deben leerse como `NH₂`/`H₂N`/`NH` horizontales.

## Alcance

- Restringir la estrategia vertical (`above`/`below`) del glifo H implícito a los
  contextos nitrogenados con contrato visual propio: `amide`, `urea` y `aromatic`
  (además de `other` para no alterar comportamiento no contratado).
- Para el contexto `amine` (definido por conectividad y órdenes de enlace en
  `_nitrogen_geometry_context`, no por proximidad geométrica) limitar los candidatos
  de `_choose_implicit_hydrogen_position` a `left`/`right`.
- Refrescar la reposición del H cuando cambia un enlace **vecino** que altera la
  clasificación del N (p. ej. dibujar `C=O` sobre el carbono ya unido al N): los N
  adyacentes a los extremos del enlace modificado se refrescan en
  `add_bond_item`/`update_bond_item`/`remove_bond_item`.
- Actualizar el contrato de amina secundaria en
  `tests/test_amide_urea_2d_continuity.py` (de `above` a horizontal) y añadir
  `tests/test_amine_label_orientation.py` con los casos etilamina, anilina,
  dietilamina, la reproducción exacta del dicetona diamino, controles de acetamida,
  urea y amonio cargado, más refrescos tras edición/undo-redo/arrastre/cambio de tema.

## Invariantes / fuera de alcance

Cambio exclusivamente de presentación. No se tocan identidad/conectividad molecular,
conteos de H, valencias ni cargas, fórmula, longitudes ni ángulos de enlace reales,
estereoquímica, persistencia `.cmsn`, undo/redo, generación 3D, RDKit, Clean2D ni
Live MM. El recorte de enlaces (`label shrink`) y el anclaje óptico del N permanecen.
No se cambia el tratamiento de amonio/N cargado, N aromático en anillo ni centros de
coordinación. No se introducen clasificadores duplicados, dependencias nuevas ni
nuevos paquetes.

## Módulos probables

`src/chemuson/gui/canvas/canvas_text.py` (política de `separated_nh` y candidatos),
`src/chemuson/gui/canvas/canvas_bond_geometry.py` (helper de refresco de contexto),
`src/chemuson/gui/canvas/canvas_structure.py` (llamadas de refresco en alta/baja/cambio
de enlace), pruebas de etiquetas. Sin dependencias ni paquetes nuevos.
