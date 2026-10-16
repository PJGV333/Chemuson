# Validación — robustez de depiction y dibujo

## Estado

Están implementados los ajustes medidos de etiqueta↔enlace, el acceso de dibujo en bordes y la geometría Haworth C3. Clean2D permanece sin cambios por falta del documento exacto que reproduce el caso del propietario. OpenSpec/baseline inicial comprometido en `632f576`; el baseline global no está completamente verde: `pytest -q` agotó 480 s y el shard agregado 4 agotó 300 s. Los 2267 node IDs baseline quedaron cubiertos en particiones terminadas: 2247 passed, 20 skipped, 0 fallos de aserción. El verificador oficial baseline de ocho shards no se declara success; detalle en `baseline.md`.

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

Baseline Qt: `_is_on_paper` con `PAPER_MARGIN=40` rechazó QTest de enlace desde x=`20` (dentro del papel) y x=`-40` (borde visible); un click no aceptado dejó grafo/undo vacíos. El `sceneRect` inicial `(-884,-884,2568,2768)` no incluía un átomo en x=`-1124`. Se retiraron ambos gates de página (herramienta y plantilla); ahora los límites incluyen coordenadas moleculares más el viewport visible y se actualizan al cambiar el índice de Undo. Un QTest arranca en x=`-40` y suelta en x=`-80`, fuera del viewport: crea el enlace y expande el límite izquierdo de escena; Undo/Redo conserva el grafo y Undo no salta el foco. QTest de anillo confirma ciclo 6/6 y grado 2 por átomo en el borde.

La sonda Haworth baseline detectó OH de C3 superpuesto con tres trazos de anillo/CH₂OH. Se reorientó sólo esa coordenada; la prueba cubre α/β × bold/plain, cero colisiones glyph↔glyph/glyph↔bond, 11 átomos y 11 enlaces invariantes. La inserción Haworth en el borde y Undo/Redo también pasan. Imágenes antes/después exploratorias en `/tmp/chemuson-depiction-campaign-baseline/haworth-chair-*.png`; owner visual review sigue pendiente.

## Ruff / evidencia

Ruff focal en los archivos de esta fase, compileall focal y `git diff --check` pasaron. El único error baseline global Ruff sigue siendo el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó. Exploraciones bajo `/tmp`; `tools/visual/render_bond_junction_gallery.py` reproduce la única lámina comparativa versionada.

## Validación final

La colección actual contiene **2275 node IDs** (8 shards: `[285, 285, 285, 284, 284, 284, 284, 284]`, SHA256 `b37eec04282f93bea22acfd35e086f110094278a8d29d476ac7bdb7ecdf73319`). `pytest --collect-only -q` y la generación del plan coincidieron con el manifiesto actualizado.

- Focales de canvas, selección, plantilla Haworth y wedges: **73 passed**; arquitectura: **280 passed**; workflow del manifiesto/shards: **4 passed**.
- `python -m compileall -q src tests tools packaging`, Ruff focal con `F401,F811,F821,E722,E741`, `git diff --check` y `openspec validate harden-molecular-depiction-and-drawing --strict`: success.
- Ruff global selectivo conserva un único error baseline, `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se cambió.
- Shards terminados: 0 `284 passed, 1 skipped`; 1 `285 passed`; 2 `281 passed, 4 skipped`; 3 `279 passed, 5 skipped`; 5 `282 passed, 2 skipped`; 6 `280 passed, 4 skipped`; 7 `284 passed`. El shard 4 excedió el límite de 300 s al 95% (sin fallo de aserción). Sus 284 IDs se recuperaron en tres particiones: Clean2D `86 passed`, resto no UI `130 passed, 4 skipped`, y `test_ui_tool_rail.py` `64 passed`.
- Cobertura total reconciliada: **2255 passed, 20 skipped, 0 fallos de aserción**. Aun así, el verificador oficial de ocho shards no se declara verde por el timeout del shard 4; el baseline monolítico `pytest -q` también agotó 480 s y no se repitió como corrida agregada.

## Entrega

**NOT READY** para cerrar la campaña: falta el documento/caso exacto del fallo Clean2D del propietario y su revisión visual manual; además, la corrida oficial de ocho shards tuvo el timeout documentado. Los cambios y validaciones quedan en `feat/molecular-depiction-excellence`; sólo commits locales, sin push no autorizado, merge a `main`/`release/v0.3.0-beta.1-prep`, tags ni releases. La única lámina comparativa versionada es `docs/visual/bond-junctions-before-after.png`; capturas exploratorias permanecen en `/tmp`.
