# Proposal: QA final y cierre de release de la modernización UI (Fase 8)

## Propósito

Las Fases 1–7 de modernización UI están integradas en `main` (`140a080`) y la revisión manual en KDE/Wayland fue aprobada. La Fase 8 demuestra que esa integración es un Release Candidate sano, comprueba el mecanismo de distribución existente, genera evidencia final, actualiza documentación y archiva formalmente los OpenSpec completados.

Esta fase es **QA + release closure exclusivamente**. No añade funcionalidad ni altera comportamiento de producto.

## Alcance

- Capturar baseline verificable de `main` y comparar la suite completa con los fallos históricos documentados.
- Ejecutar QA arquitectónica, Ruff scoped, pruebas UI dirigidas y smoke Qt light/dark, tamaños normal/compacto e HiDPI.
- Auditar el empaquetado soportado actualmente por el repositorio y probar instalación/artefacto fuera del checkout, sin introducir otra tecnología.
- Generar un conjunto pequeño de capturas before/after y README de evidencia.
- Actualizar `docs/ui-modernization/PLAN.md`, la documentación de usuario/release notes conforme a la convención del repo y declarar el resultado de la auditoría de versión.
- Archivar Fases 1–7 exclusivamente con el comando oficial OpenSpec; mantener activo el OpenSpec de esta fase hasta completar sus tareas.

## Fuera de alcance

- Cualquier nueva funcionalidad, Fase 9 o refactor oportunista.
- Cambios en Clean2D, ChemName, química/geometría de plantillas, Fischer, Haworth, silla, tetrandrina, canvas geometry, hit-testing, molecular graph, molblocks o persistencia `.cmsn`.
- Añadir un nuevo sistema/dependencia de packaging o hacer bump de versión sin contrato de release confirmado.
- Corregir fallos históricos ajenos. Si QA encuentra una regresión nueva, se detiene y se reporta.

## Impacto previsto

Principalmente OpenSpec y documentación (`openspec/changes/`, `docs/ui-modernization/`, manual y release notes existentes), más aproximadamente cinco imágenes finales en `docs/ui-modernization/after/`. Una corrección de código sólo se considerará si se demuestra pérdida de assets en el mecanismo oficial de distribución; cualquier problema ajeno a la modernización detiene el trabajo.

## Identificador

OpenSpec CLI 1.5.0 exige que el nombre de un cambio comience con una letra. Por eso se usa `ui-modernization-final-qa-release-2026-10-02` en vez del ejemplo con fecha inicial.