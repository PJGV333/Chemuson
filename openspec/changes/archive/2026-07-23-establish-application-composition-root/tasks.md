## 0. Baseline

- [x] 0.1 Capturar estado Git, `compileall`, colección, suite, Ruff y OpenSpec; aceptación: resultados exactos y limitaciones ambientales quedan en `baseline.md`.

## 1. Reconocimiento

- [x] 1.1 Inventariar consumidores de `run_app`, imports de startup, entry points de packaging y dependencias directas M08/M19; aceptación: no existe consumidor de compatibilidad que obligue a conservar `gui.main_window.run_app`.

## 2. Caracterización

- [x] 2.1 Añadir pruebas no Qt del CLI y orden de lifecycle; aceptación: `--version`/ayuda no cargan GUI y fakes registran la secuencia histórica de arranque, recuperación, event loop, errores y exit code.

## 3. Paquete de composición

- [x] 3.1 Crear `src/chemuson/app/` como parte de M19; aceptación: su `__init__.py` es liviano y no importa PyQt6 ni GUI.
- [x] 3.2 Trasladar `run_app()` al módulo interno de bootstrap; aceptación: conserva literalmente orden, metadata, autosave check, mensajes de fallo y terminación.

## 4. Entry point mínimo

- [x] 4.1 Cambiar `chemuson.__main__.main` para delegar mediante import diferido al composition root; aceptación: `chemuson`, `python -m chemuson` y `--version` conservan su contrato.
- [x] 4.2 Retirar `run_app()` de `gui/main_window.py`; aceptación: M08 sólo define la ventana y no conserva alias ni reexport de bootstrap.

## 5. Catálogo y límites

- [x] 5.1 Actualizar M19 con paths, dependencias, API y tests reales; aceptación: M19 posee `app/`, depende de M08/M15/M18 y ningún módulo M00-M18 depende de M19.
- [x] 5.2 Recalcular M08 después del movimiento; aceptación: sólo se eliminan dependencias que ya no tengan imports directos y no se añade ninguna excepción o ciclo.
- [x] 5.3 Reforzar pruebas AST de ownership, dirección e import isolation; aceptación: un import M00-M18→M19, un bootstrap dentro de `main_window.py` o una carga Qt en la ruta CLI ligera falla.

## 6. Documentación

- [x] 6.1 Crear `docs/modules/M19-bootstrap.md` y actualizar índice/arquitectura; aceptación: documentan el composition root real, sus dependencias y el límite con el futuro application shell.
- [x] 6.2 Confirmar la deuda arquitectónica en cero; aceptación: documentación y catálogo no contienen afirmaciones históricas de excepciones activas.

## 7. Validación

- [x] 7.1 Ejecutar pruebas focalizadas y arquitectónicas, `compileall`, OpenSpec estricto, suite completa disponible y Ruff; aceptación: resultados y límites del entorno quedan documentados en `validation.md`, sin ocultar fallos.

## 8. Revisión y archivo

- [x] 8.1 Revisar diff y alcance antes de archivar; aceptación: no hay cambios en menús, canvas, items, controllers, serialización, renderizado ni lógica química.
- [x] 8.2 Archivar sólo después de aprobación explícita; aceptación: specs se fusionan y OpenSpec global vuelve a validar estrictamente.
