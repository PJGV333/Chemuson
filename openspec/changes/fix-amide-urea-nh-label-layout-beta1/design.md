# Diseño — ancla N fija y glifo H independiente

## Diagnóstico visual previo

Las capturas `openspec/changes/improve-amide-urea-2d-continuity/evidence/after/urea-after-light-1x.png` y `.../urea-after-dark-1x.png` ya incluyen el arreglo anterior del trim. En ambas, el N/H sigue apareciendo como un único texto horizontal y el H se fuerza al costado inline; `_base_label_offset()` además desplaza el label completo en heteroátomos con hidrógeno. Por tanto, no se debe tocar el grafo o el algoritmo previo de trimming: el defecto restante es la colocación del glifo H con respecto al N.

## Decisiones

1. Mantener `_build_atom_label()` como fuente de `NH`/`HN`/`NH2`/`H2N` y conteo; no inferir química a partir del texto. Para render, `AtomItem` divide solo las etiquetas de N neutro no coordinado en N anclado y H/H₂ visual independiente; `toPlainText()` y la etiqueta lógica siguen exponiendo el label completo. N cargado/coordinado conserva íntegro su layout probado de cuñas/estilos.
2. Anclar el centro de tinta del glifo N al centro del `AtomItem`. En N con H inline, `_base_label_offset()` y el desplazamiento extra por longitud custom no desplazan el N. No se escriben coordenadas al modelo.
3. Seleccionar entre derecha, izquierda, arriba y abajo mediante el bisector del espacio angular libre de enlaces estructurales; puntuar adicionalmente las formas visibles de enlaces, etiquetas y cargas cercanas con un margen óptico acotado. Si todas las opciones chocan, elegir la de menor colisión. No se rotan glifos.
4. La forma/bounds del label principal incorpora el H/H₂ separado, de modo que trimming, colisiones y exports que consultan `label.shape()` siguen viendo toda la tinta. El clearance y estilo actuales permanecen sin cambios.
5. Las etiquetas N sin H y los heteroátomos O/S con H conservan su comportamiento. El layout se recalcula en `_refresh_atom_label()` después de cambios de grafo, estilo o label.
6. QTest crea las estructuras mediante herramienta de anillo/enlace del canvas; preview/commit, centros, distancias, ángulos, undo/redo y CMSN se verifican. Las imágenes son evidencia Qt offscreen, no aceptación manual nativa.

## Baseline

Estado limpio en la rama base `release/v0.3.0-beta.1-prep`, HEAD `409581dc685fd166428bde47d8355633a8eb05da`; nueva rama aislada creada después de confirmar ausencia de cambios locales: `fix/amide-urea-nh-label-layout-beta1`. Salidas exactas del baseline obligatorio: `baseline/01-git-status.log` a `baseline/05-ruff.log`; los logs grandes de colección/pruebas se comprimen sin alterar su contenido.

## Riesgos

- Cambiar solo la pintura y no `shape()` haría que el recorte de enlace ignore el nuevo H; la prueba de intersección cubre este contrato.
- Desplazar el label completo para evitar colisión incumpliría el anclaje N; se verifica con el centro de tinta frente al centro químico.
- El margen geométrico puede preferir una orientación menos natural aun libre; las capturas referencia y controles de aminas guían el ranking.
- El hidrógeno `2` mantiene el formato subíndice vigente y no se interpreta como elemento independiente.
