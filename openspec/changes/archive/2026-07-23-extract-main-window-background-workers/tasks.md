## 0. Baseline

- [x] 0.1 Capturar Git, compileall, pytest, Ruff y OpenSpec; aceptación: límites del entorno documentados.

## 1. Reconocimiento y caracterización

- [x] 1.1 Confirmar que sólo `main_window.py` consume ambos workers y que sus backends se importan dentro de `run()`.
- [x] 1.2 Añadir pruebas AST de ownership, señales, imports diferidos y ausencia de ciclo.

## 2. Extracción

- [x] 2.1 Crear `gui/background_workers.py` dentro de M08.
- [x] 2.2 Trasladar literalmente `_DescriptorWorker` y `_NameToStructureWorker`.
- [x] 2.3 Mantener thread lifecycle, jobs, cancelación y callbacks en `ChemusonWindow`.

## 3. Arquitectura y documentación

- [x] 3.1 Registrar el módulo interno y su test en M08.
- [x] 3.2 Actualizar la ficha M08 sin crear API pública.

## 4. Validación

- [x] 4.1 Ejecutar caracterización, arquitectura, compileall, Ruff, OpenSpec y suite disponible.
- [x] 4.2 Confirmar que el diff no altera UI, canvas, química, persistencia ni formatos.

## 5. Revisión y archivo

- [x] 5.1 Publicar para prueba manual; archivar sólo después de aprobación explícita.
