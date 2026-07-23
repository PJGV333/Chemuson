## 0. Baseline

- [x] 0.1 Capturar estado Git, `compileall`, colección, Ruff y OpenSpec; aceptación: resultados exactos y límites ambientales quedan documentados.

## 1. Reconocimiento

- [x] 1.1 Inventariar regiones, orden, colaboradores y atributos creados por `ChemusonWindow.__init__`; aceptación: el límite distingue composición de handlers.

## 2. Caracterización

- [x] 2.1 Añadir pruebas AST del seam y ownership; aceptación: una composición duplicada, dependencia circular o ruta fuera de M08 falla.

## 3. Application shell

- [x] 3.1 Crear `src/chemuson/gui/shell/` como implementación interna de M08.
- [x] 3.2 Trasladar literalmente estado, controllers, tabs, docks, toolbars, status bar y wiring del constructor.
- [x] 3.3 Reducir `ChemusonWindow.__init__` a inicialización base y delegación única.

## 4. Catálogo y documentación

- [x] 4.1 Actualizar M08 con ownership interno y pruebas reales.
- [x] 4.2 Documentar el límite entre composition root, application shell, ventana y canvas.

## 5. Validación

- [x] 5.1 Ejecutar caracterización, pruebas arquitectónicas, `compileall`, Ruff, OpenSpec y suite disponible.
- [x] 5.2 Revisar que no haya cambios visuales, de handlers, canvas, química, persistencia o formatos.

## 6. Revisión y archivo

- [x] 6.1 Archivar sólo después de prueba manual y aprobación explícita.
