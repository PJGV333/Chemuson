## 0. Baseline

- [x] 0.1 Capturar Git, compileall, pytest, Ruff y OpenSpec; aceptación: límites del entorno documentados.

## 1. Reconocimiento y caracterización

- [x] 1.1 Inventariar implementaciones y consumidores de las siete reglas.
- [x] 1.2 Añadir pruebas AST de ownership, aliases, wrapper e imports permitidos.
- [x] 1.3 Añadir regresiones numéricas directas para fórmulas, tolerancias y límites.

## 2. Extracción

- [x] 2.1 Crear `gui/canvas/selection_geometry.py` dentro de M09.
- [x] 2.2 Trasladar literalmente las siete reglas deterministas.
- [x] 2.3 Conservar seis aliases privados y el wrapper de grosor en `CanvasSelectionMixin`.

## 3. Arquitectura y documentación

- [x] 3.1 Registrar el módulo interno y sus pruebas en M09.
- [x] 3.2 Actualizar la ficha M09 sin crear API pública ni dependencias nuevas.

## 4. Validación

- [x] 4.1 Ejecutar regresiones focalizadas, arquitectura, compileall, Ruff y OpenSpec.
- [x] 4.2 Ejecutar colección y suite globales y comparar con el baseline EGL.
- [x] 4.3 Confirmar que el diff no altera eventos, comandos, items, química, renderizado ni persistencia.

## 5. Revisión y archivo

- [ ] 5.1 Publicar para prueba manual; archivar sólo después de aprobación explícita.
