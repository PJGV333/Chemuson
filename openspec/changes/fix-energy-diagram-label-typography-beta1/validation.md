# Validación — tipografía de menús de diagramas de energía

## Cambio

- Se añadieron helpers de presentación específicos para etiquetas de menú de presets semánticos y presets de cajas de energía.
- `QAction` sigue capturando/emitiendo `preset.name`; `energy_diagram_display_name`, los `display_name`/`title` de los presets, `EnergyDiagramPreset.label`, `kind`, y todos los parámetros de construcción permanecen intactos.
- La paleta muestra `sp²`/`sp³` en sus glifos y `Híbrido sp`, `Híbrido sp²`, `Híbrido sp³` en tooltips. Los submenús muestran las fórmulas y configuraciones d con Unicode. El formato conserva `2s`/`2p` en línea base.
- No se modificaron builders científicos, ocupaciones, electrones, energía, degeneración, orden, geometría, persistencia ni dependencias.

## Pruebas y verificaciones

- Baseline focal anterior al cambio: **86 passed** (log `baseline/06-focused-tests.log`).
- Verificación final focal de builders, presets, persistencia semántica, menús/flyouts y shell compacto: **126 passed in 218.53s**.
- Incluye selección Qt real del flyout para `hybrid_sp3`; selección de la QAction `O₂²⁻` emite exactamente `mo:O2^2-`; inventarios e IDs originales exactos; etiquetas completas de los presets MO/campo ligando; temas claro/oscuro y métricas Qt para bounds/soporte de glifos.
- SHA-256 canónico de `to_json_dict()` de los 14 presets MO y campo ligando coincide con los valores capturados sobre la beta base antes del cambio. La evidencia `semantic-preset-fingerprints.json` registra además la baseline de los 11 presets atómicos. Los títulos y labels persistentes siguen siendo ASCII de baseline, incluyendo `O2^2- MO Diagram`, `sp3` y `2s/2p`.
- `python -m compileall src tests tools packaging` — PASS.
- Ruff scoped a los seis archivos Python de producto/prueba de esta campaña (`F401,F811,F821,E722,E741`) — PASS.
- Ruff de repositorio completo conserva el único error basal no relacionado: F401 `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; el archivo no se modificó.
- `tests/test_test_workflow_contract.py` — 4 passed; node-ID manifest contract permanece válido sin añadir/quitar tests.
- `openspec validate fix-energy-diagram-label-typography-beta1 --strict` — PASS.
- `git diff --check` — PASS.
- La suite monolítica no se repitió según la autorización de alcance; su timeout y los fallos Qt/SIP conocidos están registrados en el OpenSpec de preparación beta y referenciados en `baseline.md`.

## Evidencia visual

- `evidence/before/`, `evidence/after/` y `evidence/comparison/` contienen las paletas y submenús MO/campo ligando, antes/después, claro/oscuro y comparativas lado a lado. Las capturas a escala Qt 1.5× también están en los subdirectorios `scale-1.5/`.
- Captura desde la sesión KDE disponible usando Qt 6/QPA `xcb` sobre `HDMI-A-1`; DPR 1.0 y Qt `QT_SCALE_FACTOR=1.5` (DPR 1.5, área lógica disponible 1067×553). No se cambió la configuración del escritorio.
- Qt reportó la cobertura de `₂`, `⁻`, `²`, `³`, `⁴`, `⁵`, `⁶` y `⁸`. Todas las acciones capturadas tienen al menos 20 px de margen de texto. El menú MO estrecho mide 88 px y el de campo de ligandos 199 px en escala 1×; todos caben sin recorte en las capturas a 1.5×.
- Las imágenes son capturas de widgets QMenu reales renderizados por Qt y evidencia de esta sesión; no sustituyen una aceptación visual manual adicional del propietario.
