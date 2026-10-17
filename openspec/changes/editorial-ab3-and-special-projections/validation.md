# Validación — Preset editorial AB3 y proyecciones especiales

## Alcance y decisiones

- OpenSpec creado y validado antes de implementar en `b572d9d`; baseline previo al cambio: `70efde4`. Implementación: `b1ec626`; aclaración de cohorte AB3: `27f39ce`.
- El preset `ChemUSON Editorial AB3` contiene exactamente los valores aprobados, se activa de manera opt-in desde Working/UI o mediante `--preset`, y mantiene intactos Current y los defaults de producción. La cohorte editorial es explícitamente **A1, A2, A3, B1, B2 y B3**. A4 y B4 quedan fuera; B4 continúa renderizable como fixture geométrico no editorial y está excluido de selección/calibración.
- La geometría pura y tipada está separada en `src/chemuson/gui/projections/` y registrada como subpaquete interno de M08. Haworth usa parámetros explícitos para familia/tamaño de anillo, proporción y sesgo, anómero, orientación D/L, grupos/sustituyentes y aristas frontales; Fischer y silla tienen constructores independientes. `templates.py` conserva las entradas históricas como adaptadores. La geometría no introduce `BondStereo` ni estereoquímica química implícita.
- No se cambiaron los defaults globales de dibujo, Clean2D, ChemName, ChemIO, formato `.cmsn` ni dependencias.

## Archivos incluidos

- Arquitectura/código: `architecture/modules.yml`; `src/chemuson/gui/templates.py`; `src/chemuson/gui/projections/{__init__,common,haworth,fischer,chair}.py`.
- Laboratorio/exportación: `tools/visual/depiction_lab.py`; `tools/visual/depiction_lab_app/{README.md,exporting.py,named_presets.py,profiles.py,scenes.py,ui.py}`.
- Pruebas y manifiesto: `tests/architecture/test_special_projection_boundary.py`; `tests/ci/expected_pytest_nodeids.txt`; `tests/probes/depiction_lab_ui.py`; `tests/test_depiction_lab_export.py`; `tests/test_depiction_lab_renderer.py`; `tests/test_editorial_depiction_preset.py`; `tests/test_special_projection_geometry.py`; `tests/test_wedge_double_label_junctions.py`.
- OpenSpec: `.openspec.yaml`, `baseline.md`, `baseline-command-output.txt`, `proposal.md`, `design.md`, `tasks.md`, `validation.md` y las especificaciones `molecular-depiction-visual-lab/spec.md` y `special-projection-geometry/spec.md`.
- Evidencia añadida bajo `evidence/`: comparación Haworth anterior/nueva y bundles con PNG/SVG, perfil y reporte para ambos estados, galería AB3 y exportación C3 con AB3.

## Verificación

- Pruebas focales tras aclarar la cohorte: `QT_QPA_PLATFORM=offscreen pytest -q tests/test_depiction_lab_renderer.py tests/test_depiction_lab_export.py tests/test_editorial_depiction_preset.py tests/test_depiction_lab_ui.py` — **39 passed**.
- Geometría/proyección, frontera arquitectónica, junctions Haworth y UI/exportación: **18 passed** en ejecución focal previa; el conjunto completo también está incluido en los shards finales.
- Arquitectura: `QT_QPA_PLATFORM=offscreen pytest -q tests/architecture` — **283 passed**.
- `python -m compileall src tests tools packaging` — **PASS**.
- Ruff dirigido a todos los archivos Python afectados con `F401,F811,F821,E722,E741` — **PASS**. El comando global requerido conserva sólo el F401 preexistente `import math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`, registrado en baseline y fuera de alcance.
- Plan final: **2,341 node IDs**, reconciliados con `tests/ci/expected_pytest_nodeids.txt`; los ocho shards se ejecutaron con límites de 300 s y `tools/ci_pytest_shards.py verify` confirmó cobertura única/completa. Resultado agregado: **2,321 passed, 20 skipped, 0 failed**.
- `openspec validate editorial-ab3-and-special-projections --type change --strict --no-interactive` — **PASS**.
- `git diff --check` — **PASS**.
- Baseline monolítico: `pytest -q` había registrado un fallo Qt conocido y dependiente del orden en `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` (más 2,305 passed y 20 skipped). No se repitió la suite monolítica; la cobertura final completa se verificó por shards, donde ese test pasó. Se conserva aquí el fallo basal sin ocultarlo.

## Evidencia visual reproducible

Los exports se generaron en Linux con Qt/PyQt 6.11.0, tema claro y **working tree limpio**, desde renderer SHA `27f39ce669bf794db1bad96c0a404c27885bb0ff`.

- [Comparación antes/después de Haworth](evidence/haworth-before-after.png) — SHA-256 `539581f7b9bf84d286e0cceb977745d4d688465b6decd059be97bb0895be0030`. Compara el PNG de baseline `70efde4` con el nuevo C3 usando Current, sin cambios en parámetros de estilo.
- [PNG Haworth anterior](evidence/haworth-before/current-C3-light.png) — `58b59de56bc1e77d5230d8f942cd9e5fb3559ecc739162fe18b1a3a7c84a0c48`; [PNG nuevo](evidence/haworth-after/current-C3-light.png) — `ba9782f1641713a9496d4aa42b75a9120f57da920704c866f77be348962fe0ad`. Reports y perfiles acompañantes registran procedencia y parámetros.
- [Galería completa Current/AB3, 12 casos](evidence/editorial-ab3-gallery/contact-sheet-light.png) — `9112cd481155b187a293abfbc2f65a04ad727cbb67c7b7806b019e123ca519d2`. A4 está marcado diagnóstico/excluido; B4 está rotulado fixture geométrico/no editorial/excluido. El reporte enumera los seis casos elegibles. El perfil reproducible está junto a la imagen.
- [Comparación de exportación C3 con AB3](evidence/editorial-ab3-haworth/comparison-C3-light.png) — `7ecf1012f47daf5dad5c8386819440fa58e063f2125109cef5328725570df776`; incluye reporte y perfil. La exportación se reproduce con `python tools/visual/depiction_lab.py --export-case C3 --output <dir> --preset "ChemUSON Editorial AB3" --theme light`.

La evidencia permite revisión visual; no promueve el preset a defaults ni cambia la aprobación de valores owner-approved.

## Ruff y límites

El Ruff global exacto falla únicamente por el F401 basal documentado; no se modificó ese test ajeno al alcance. La ejecución `pytest -q` monolítica basal también conserva el fallo Qt ya descrito. La cobertura shard final sí quedó completa y verde. No se afirma promoción a producción ni se alteró el estado químico de los grafos de proyección.
