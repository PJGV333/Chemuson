# Validación — robustez de depiction y dibujo

## Estado

La corrección de trimming etiqueta↔enlace está implementada; Clean2D y el input del canvas siguen sin cambios de producción. OpenSpec/baseline inicial comprometido en `632f576`; OpenSpec estricto pasó. El baseline global no está completamente verde: `pytest -q` agotó 480 s y el shard agregado 4 agotó 300 s. Los 2267 node IDs quedaron cubiertos en particiones terminadas: 2247 passed, 20 skipped, 0 fallos de aserción. El verificador oficial de ocho shards no se declara success; detalle en `baseline.md`.

## Diagnóstico Clean2D (sin cambio)

Sondas Qt del motor y del controlador, con `QT_QPA_PLATFORM=offscreen`, registradas bajo `/tmp/chemuson-depiction-campaign-baseline/`:

- Ciclo manual cerrado de seis carbonos (simple no aromático, Kekulé no aromático y aromaticidad explícita): quick/publication conservaron la topología; casos distorsionados seleccionaron `rdkit_isolated` o `simple_aromatic_template`, mantuvieron cero cruces y dieron razón de degeneración `0.9069` (antes `0.7528`). La snapshot de grafo permanece idéntica.
- Una cadena visualmente casi cerrada con seis átomos pero **cinco enlaces** (terminales separados `5.83`): Ctrl+K la reacomodó como cadena/zigzag (`rdkit_isolated`), sin cambiar ninguno de sus cinco enlaces. El dibujo inicial no codifica un ciclo; este resultado no es evidencia de que Clean2D haya borrado un enlace de anillo.
- La herramienta de enlace Qt creó un ciclo manual de seis átomos/seis enlaces con grado 2 en cada átomo al cerrar sobre el átomo inicial. No se observó fallo del snap en esa interacción reproducible.
- Diagnóstico no concluido para el documento reportado por el propietario: se desconoce si su grafo tenía el sexto enlace. No se va a inferir ni agregar conectividad desde Clean2D. Mantener el caso exacto pendiente hasta reproducir un ciclo de seis enlaces que se colapse o conocer el archivo fuente.

## Wedge-to-label (medición base)

La galería ampliada usa la ruta de producción y cubre F/Cl/Br/N⁺/CH₂OH, extremos estrecho/ancho, sólida/hash, controles simple/doble/triple/bold, bifurcaciones y dos temas. El before baseline está en `/tmp/chemuson-depiction-campaign-baseline/label-gallery-before/`. Distancia axial de formas a 1×: baseline de cuña sólida `3.50` unidades en F/Cl/Br (hasta `9.90` en base ancha); hash `7.58–7.82` en F/Cl/Br (base ancha hasta `7.83`) y hasta `9.74` en P; N⁺ sólido `4.00`/hash `8.72`.

La implementación sustituye el bounding-ellipse genérico por soporte direccional de los glifos reales (incluida carga) y clearance de trazo; elimina el margen extra basado en ancho completo en la base sólida. Los marks hash mantienen su cantidad y espaciado uniforme; el trim se calcula para que el primer/último mark conserve clearance sin el antiguo medio intervalo vacío. La matriz Qt cubre 5 labels × 4 orientaciones × 6 estilos × ambos extremos; todos los casos no intersectan y la distancia medida está entre `0.5` y `5.0` unidades. Pruebas focales de wedge/galleries: 41 passed. Comparación expandida versionada en `docs/visual/bond-junctions-before-after.png`; capturas intermedias permanecen en `/tmp/chemuson-depiction-campaign-baseline/`. La revisión visual del propietario sigue pendiente.

## Canvas y Haworth

Causa preliminar del límite visible confirmada en código: `CanvasSelectionInputMixin.mousePressEvent` descarta herramientas de dibujo si `_is_on_paper(scene_pos)` es falso; la zona válida excluye `PAPER_MARGIN=40` px. Además, `_update_scene_rect` usa la hoja + un margen ligado al viewport, no las coordenadas del grafo. Falta la reproducción Qt de borde/Undo y el diseño de expansión de navegación.

Pendiente revisar visualmente la plantilla Haworth y geometría de cuñas en bifurcaciones antes de cambiar producción.

## Ruff / evidencia

Ruff focal en los archivos de esta fase, compileall focal y `git diff --check` pasaron. El único error baseline global Ruff sigue siendo el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó. Exploraciones bajo `/tmp`; `tools/visual/render_bond_junction_gallery.py` reproduce la única lámina comparativa versionada.

## Entrega

Pendiente. La entrega permanecerá exclusivamente en `feat/molecular-depiction-excellence`; sin merge a `main`/`release/v0.3.0-beta.1-prep`, tags ni releases. La revisión visual del propietario es independiente de las pruebas automáticas.
