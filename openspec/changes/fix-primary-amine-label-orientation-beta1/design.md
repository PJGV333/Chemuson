# Diseño — aminas ordinarias con NH₂/H₂N/HN horizontales

## Causa raíz confirmada

1. `canvas_text._refresh_atom_label` activa `separated_nh` para **cualquier** N neutro
   no centros de coordinación con etiqueta `NH/HN/NH2/H2N`, sin consultar el tipo
   químico del nitrógeno.
2. `_choose_implicit_hydrogen_position` clasifica candidatos
   `right/left/above/below`; gana el que menos colisiona y, entre empates, el más
   cercano al sector abierto. Con enlaces inclinados ±75° (repro del dicetona
   diamino) los candidatos verticales ganan: en la reproducción oficial los dos
   amino quedaron `above` y `below` respectivamente.
3. Los contratos previos de amidas/ureas (`fix-amide-urea-nh-label-layout-beta1`)
   **sí** requieren arriba/abajo en vertices conjugados; por eso no se puede
   desactivar la estrategia global.

El clasificador `_nitrogen_geometry_context` (`canvas_bond_geometry.py`) ya existía
para la selección de geometría sp2 y clasifica por conectividad y órdenes de enlace:
`coordination`, `amide`, `urea`, `aromatic`, `other`, `amine`. Anilina (N unido por
enlace simple a carbono aromático) cae correctamente en `amine`.

## Decisión

- Reutilizar el clasificador existente; **no** crear uno nuevo.
- En `_refresh_atom_label`, pasar `allow_vertical = context != "amine"` a
  `_choose_implicit_hydrogen_position`. `amide`/`urea`/`aromatic`/`other` conservan
  las cuatro direcciones exactas de antes; `amine` solo `left`/`right` (la clasificación
  depende de conectividad y órdenes de enlace, nunca de proximidad geométrica).
- Confinar el cambio a la selección de candidatos; el layout de glifos, el anclaje
  óptico y el recorte de enlaces quedan intactos.
- Consecuencia de contrato documentada: el caso parametrizado
  `secondary-amine` de `test_amide_urea_2d_continuity.py` esperaba `above`; pasa a
  `left` (H₂N→HN prefijo con H a la izquierda, horizontal). Es exactamente el
  comportamiento que la solicitud pide corregir.

## Efecto secundario detectado y corregido

La clasificación por contexto evalúa los **vecinos** del N. Antes, los refrescos de
label solo cubrían los extremos del enlace editado: al dibujar `C=O` sobre el carbono
vecino de un N, el N quedaba con posición calculada bajo el contexto obsoleto. Se
añade `_nitrogen_context_dependent_refresh_ids(endpoint_a, endpoint_b)` en
`canvas_bond_geometry.py` (N adyacentes a los extremos del enlace modificado) y se
invoca desde `add_bond_item`, `update_bond_item` y `remove_bond_item`. Sin refresco,
el contrato de acetanilida/urea dibujados fuera de orden regresaba.

## Archivos afectados

- `src/chemuson/gui/canvas/canvas_text.py`
- `src/chemuson/gui/canvas/canvas_bond_geometry.py`
- `src/chemuson/gui/canvas/canvas_structure.py`
- `tests/test_amide_urea_2d_continuity.py` (un valor de contrato)
- `tests/test_amine_label_orientation.py` (nuevo)

## Riesgos de regresión y mitigación

- Ureas/amidas con H vertical: protegidos por los contratos existentes
  (`test_amide_urea_2d_continuity.py`, `test_functional_group_hydrogen_adjustment.py`),
  todos en verde.
- N aromático (pirrol) y centros de coordinación: no se tocan; `aromatic` conserva
  cuatro direcciones y coordinación sigue excluida de `separated_nh`.
- Rendimiento: los refrescos extra iteran los enlaces del modelo al crear/eliminar
  enlaces (O(|bondos|), igual que `_refresh_bond_ring_flags` vecino); sin impacto
  medible en la batería de enlaces.
- Desplazamientos y trims: `refresh_atom_labels` solo recalcula la posición del glifo H;
  verificado geométricamente (rectángulos, distancias, solapes) en el nuevo test.

## Pruebas

Suite dirigida completa de suites de etiquetas/H/enlaces (≈249 tests en subconjuntos
`grep`-dirigidos), diagnóstico offscreen con posiciones reales de glifos y capturas
claro/oscuro en `evidence/`. Ver `validation.md`.

## Revisión 2 — secundarias en cadena y percepción de longitud

**Regresión de la regla inicial.** Restringir todo el contexto `amine` a
`left`/`right` colocaba el glifo `H` de una amina secundaria `-NH-` de cadena
**sobre el eje de la cadena**, dentro del sector de uno de sus dos enlaces. El
recorte por etiqueta (`_update_bond_label_shrink`) usa la forma visible completa
del label, que incluye el glifo H; al extenderse el H hacia el enlace, el trim de
ese extremo crecía de ~8.5 px a ~17.9 px y el trazo dibujado pasaba de 31.5 px a
22.1 px en un zigzag simétrico de enlaces de 40 unidades. El otro N–C conservaba
31.5 px: la cadena se veía desproporcionada junto a la etiqueta. Clean2D (Ctrl+K)
reubica el N y enmascara parcialmente el efecto, pero el dibujo normal nacía mal.

**Heurística final (localizada, misma función).** En `_refresh_atom_label` el
flag `allow_vertical_h` del contexto `amine` depende ahora de la etiqueta:
`NH2`/`H2N` (una sola conexión, el apilado vertical del subíndice es el defecto
original) se mantienen horizontales; `NH`/`HN` (dos conexiones, la geometría
local define el sector libre) recuperan los cuatro candidatos. La jerarquía de
`_choose_implicit_hydrogen_position` ya prioriza correctamente para secundarias:
0 colisiones primero y, entre empates, menor `direction_cost` respecto al sector
abierto `-Σ u(enlaces)`, lo que sitúa el H fuera del abanico de enlaces y deja
simétricos los trims de ambos N–C. No se añadieron puntuaciones nuevas: la
regla anterior era el único cambio y resultó excesiva.

**Contrato revocado.** El valor `left` del parámetro `secondary-amine` de
`test_qtest_amine_nh_glyphs_reposition_without_moving_the_nitrogen` vuelve a
`above` (sector abierto del V dibujado por la herramienta). Los nuevos tests de
cadena miden la longitud dibujada real de los dos N–C adyacentes (path del item
de enlace) y exigen igualdad ±0.5 px.

**Manifiesto CI.** Los node IDs de `tests/test_amine_label_orientation.py` (14)
quedaron sin registrar en la revisión 1; ahora `update-expected --confirm` lleva
el manifiesto de 2610 a 2624 y `tests/ci/pytest_file_durations.json` recibe la
entrada del archivo (`local_mean`, 0.05 s/test, conservador frente a 0.37 s/14
tests medidos en local). El contador fijado de archivos del perfil en
`tests/test_ci_pytest_shards.py` pasa de 253 a 254 y `plan` vuelve a salir 0 con
8 shards.
