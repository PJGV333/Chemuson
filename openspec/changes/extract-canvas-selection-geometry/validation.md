# Validación

Fecha: 2026-07-23

## Resultado

- Arquitectura completa más regresiones numéricas nuevas: **229 passed**.
- Regresiones nuevas: delta angular, rotación, escalado, igualdad tolerante y
  normalización de escala/grosor correctas.
- Contratos AST de la extracción: **4/4**.
- `python -m compileall -q src tests tools packaging`: correcto.
- Ruff sobre archivos tocados: correcto.
- Ruff global: sólo conserva el F401 histórico fuera del cambio.
- `openspec validate extract-canvas-selection-geometry --strict`: válido.
- `openspec validate --all --strict`: 23/23 válidos.
- `git diff --check`: correcto.

## Suite completa

- `pytest --collect-only -q`: 1070 pruebas detectadas y 56 errores conocidos
  por ausencia ambiental de `libEGL.so.1`.
- `pytest -q`: reproduce exactamente los mismos 56 bloqueos de colección.
- No apareció ningún fallo distinto del bloqueo Qt/EGL ya documentado.

## Evidencia arquitectónica

- Las siete reglas deterministas viven en
  `gui/canvas/selection_geometry.py`.
- El módulo sólo importa `math`, `__future__` y `QPointF` desde QtCore.
- `CanvasSelectionMixin` conserva seis aliases estáticos y un wrapper mínimo.
- M09 conserva sus dependencias, API pública, cero excepciones y cero ciclos.
- `canvas_selection.py` bajó de 5154 a 5128 líneas.
- La ficha M09 registra exactamente los 21 archivos Python actuales del
  paquete, incluido `__init__.py`.

## Alcance revisado

No se modificaron eventos, jerarquía de mixins, handlers, selección, drag,
hit testing, overlays, comandos undo/redo, transformación de items, ramas,
trackball 3D, química, renderizado, persistencia, formatos ni estilos.

El cambio requiere una prueba manual en Qt real antes de archivarse o
fusionarse.
