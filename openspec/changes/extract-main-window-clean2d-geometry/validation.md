# Validación

Fecha: 2026-07-23

## Resultado

- Pruebas focalizadas iniciales: **101 passed**.
- Suite arquitectónica completa más regresiones puras nuevas: **215 passed**.
- Regresiones puras nuevas: escala, centro, rotación, reflexión, traslación e
  hidrógenos faltantes correctas.
- Contratos AST de la extracción: **4/4**.
- `python -m compileall -q src tests tools packaging`: correcto.
- Ruff sobre archivos tocados: correcto.
- Ruff global: sólo conserva el F401 histórico fuera del cambio.
- `openspec validate extract-main-window-clean2d-geometry --strict`: válido.
- `openspec validate --all --strict`: 22/22 válidos.
- `git diff --check`: correcto.

## Suite completa

- `pytest --collect-only -q`: 1050 pruebas detectadas y 56 errores conocidos
  por ausencia ambiental de `libEGL.so.1`.
- `pytest -q`: reproduce exactamente los mismos 56 bloqueos de colección.
- No apareció ningún fallo distinto del bloqueo Qt/EGL ya documentado.

## Evidencia arquitectónica

- Las cinco implementaciones viven en `gui/clean2d_geometry.py`.
- El módulo sólo importa `math` y `__future__`.
- `main_window.py` conserva cinco aliases estáticos privados.
- M08 conserva sus dependencias, API pública, cero excepciones y cero ciclos.
- `main_window.py` bajó de 2319 a 2173 líneas.

## Alcance revisado

No se modificaron algoritmos, thresholds, controllers, canvas, renderizado,
química, persistencia, formatos, acciones ni estilos. El cambio requiere una
prueba manual en Qt real antes de archivarse o fusionarse.
