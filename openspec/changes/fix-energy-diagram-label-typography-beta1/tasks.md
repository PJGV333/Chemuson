# Tareas — tipografía de diagramas de energía

## A. Alcance y baseline
- [x] A.1 Verificar beta base `582fd0e2f80131a2fc9ae4e5a1a187c565bcb470`, worktree autorizado limpio y nueva rama independiente.
- [x] A.2 Inspeccionar todos los textos visibles de la paleta/submenús, tooltips y títulos persistidos relacionados.
- [x] A.3 Registrar baseline de status, compileall, collect-only, Ruff y pruebas focalizadas; referenciar los fallos monolíticos de Qt/SIP ya documentados sin reintentar la suite completa.
- [x] A.4 Capturar menú antes del cambio en tema claro/oscuro con fuente y escala registradas.

## B. Presentación
- [x] B.1 Añadir formato UI explícito para fórmulas MO, cargas, hibridación y conteos d; no inferir subíndices/superíndices por regex global.
- [x] B.2 Usar etiquetas Unicode visibles en las acciones/tooltips de los menús y paleta correspondiente.
- [x] B.3 Conservar IDs, objetos/propiedades de preset, texto de diagrama, títulos y datos persistidos.

## C. Regresiones y evidencia
- [x] C.1 Probar el catálogo original completo, las etiquetas exactas y las señales de selección (`mo:O2^2-`, `hybrid_sp3`).
- [x] C.2 Comparar datos de los diagramas construidos antes/después para B₂/C₂/N₂/O₂/F₂ y los restantes presets relevantes.
- [x] C.3 Probar persistencia y que `2s`/`2p` quedan inalterados; inspeccionar el render/bounds de glifos sin recorte.
- [x] C.4 Capturar comparaciones de submenús en temas claro/oscuro y escalas disponibles; anotar el límite de aceptación visual nativa si aplica.

## D. Verificación y entrega
- [x] D.1 Ejecutar tests focalizados, `compileall`, Ruff scoped, validación OpenSpec y `git diff --check`.
- [x] D.2 Revisar diff y artefactos, crear commit descriptivo y dejar el worktree limpio.
- [ ] D.3 Publicar solo esta rama con push normal si hay autenticación; verificar igualdad de SHA remoto/local, sin merge, PR, tag ni Release.
