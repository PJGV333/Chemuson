# Baseline

Fecha: 2026-07-23

- Rama: `architecture/phase10-extract-canvas-selection-geometry`.
- Base remota integrada: `main` `a26fa638`.
- Cierre formal de fase 9 aislado en `6cd6c66`.
- `git status --short`: limpio antes de crear los artefactos de fase 10.
- `python -m compileall -q src tests tools packaging`: correcto.
- `pytest --collect-only -q`: 1050 pruebas recolectadas y 56 errores de
  importación por ausencia ambiental de `libEGL.so.1`.
- `pytest -q`: reproduce los mismos 56 errores durante colección.
- Ruff global: reproduce únicamente el F401 histórico de
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`.
- `openspec validate --all --strict`: 22/22 válidos tras archivar la fase 9.

## Reconocimiento

- `canvas_selection.py`: 5154 líneas antes de la extracción.
- Cinco reglas candidatas ya son `staticmethod`.
- `_normalize_label_scale` no usa estado de instancia.
- `_normalize_custom_stroke` sólo lee el grosor heredado actual.
- `_optional_float_equal` y `_point_equal` también son consumidas por otros
  mixins a través del MRO, por lo que se preservarán sus nombres privados.
