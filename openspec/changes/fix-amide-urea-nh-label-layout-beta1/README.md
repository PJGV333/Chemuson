# Disposición visual N–H en ChemUSON beta.1

- Alcance y contexto: `proposal.md`, `design.md`, `tasks.md`, `specs/nh-label-layout/spec.md`.
- Baseline exacto: `baseline.md` y `baseline/`.
- Capturas QTest anteriores y métricas: `evidence/before/`.
- Capturas posteriores y métricas: `evidence/after/`.
- Comparación lado a lado en los temas claro/oscuro: `evidence/comparison/`.
- Reproducción: `evidence/before/capture_nh_layout_before.py` se ejecuta contra el checkout base; `evidence/after/capture_nh_label_layout_after.py` contra el código actual. `evidence/compare_layout_captures.py` compone las comparaciones.

Las capturas son QImage/Qt offscreen. La página del documento permanece blanca en ambos temas; el layout del glifo y los contrastes N/H se verifican en ambos casos. No equivalen a aceptación visual manual en un escritorio nativo.
