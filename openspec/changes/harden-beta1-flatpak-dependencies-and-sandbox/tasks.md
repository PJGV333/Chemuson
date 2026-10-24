# Tareas — dependencias y sandbox de beta.1

## A. Auditoría y baseline
- [x] A.1 Confirmar HEAD/árbol limpio, remotos y nueva rama aislada desde `08a90c199f526639453a8ae8b82da957a1458d08`.
- [x] A.2 Capturar compileall, colecta, Ruff y status previos a cambios; no ejecutar la suite completa.
- [x] A.3 Auditar imports/usos de Pillow, APIs Qt de imagen y rutas Flatpak; consultar OSV para Pillow 12.0.0, 12.3.0 y 13.0.0.
- [x] A.4 Comprobar QImage, `rdMolDraw2D` y tests de persistencia de imágenes con PIL bloqueado.

## B. Cambio acotado
- [x] B.1 Eliminar Pillow como dependencia directa/módulo Flatpak, excluir recursos PIL transitivos de PyInstaller y añadir validación CArchive; corregir documentación/tests.
- [x] B.2 Sustituir el permiso global de home por permisos de datos/configuración ChemUSON; conservar y explicar el permiso de red.
- [x] B.3 Añadir contrato Flatpak/dependencias sin añadir comportamiento químico ni cambiar la suite completa.
- [x] B.4 Registrar riesgos residuales y estados de aceptación exacta como pendientes; conservar los gates de publicación.

## C. Verificación
- [x] C.1 Ejecutar tests focalizados de dependencia, imagen, manifiesto, distribución, ChemName/RDKit y smoke 2D→3D apropiado.
- [x] C.2 Ejecutar compileall, OpenSpec estricto, Ruff focal y `git diff --check`; no disparar `test.yml`.
- [x] C.3 Validar el manifiesto Flatpak y revisar diff/estado; commit/push normal sólo en la rama aislada.
- [x] C.4 Ejecutar Build Preview de cuatro formatos desde SHA `e4710b938f0096ff2d02f35010117bb9a6ed8286`; comprobar jobs, `publication=false`, hashes, provenance, ausencia de recursos top-level PIL en CArchive/Flatpak, RDKit/ChemName y smoke 2D→3D.
- [x] C.5 Anotar resultados exactos y artefactos en OpenSpec/notas/matriz sin convertir smoke en aceptación manual.

## D. Dictamen histórico al cerrar el cambio de seguridad (2026-10-23)
- [x] D.1 Verificar que `main` y release-prep no se modificaron, y que no había tag/Release/publicación en ese cierre.
- [x] D.2 Dictaminar **NO-GO en ese momento**; no crear tag mientras faltara aceptación/validación interactiva o controles de publicación.

## E. Integración autorizada posterior (2026-10-24)
- [x] E.1 Integrar los cuatro commits de esta rama en `release/v0.3.0-beta.1-prep` mediante fast-forward-only; mantener preservados los hashes originales.
- [x] E.2 Registrar la aceptación manual limitada y la autorización posterior del propietario sin reasignar PASS a funciones no verificadas ni a artefactos de SHA desconocido.
- [ ] E.3 Verificar Build Preview y paquete oficial generados desde el SHA definitivo; no reutilizar #39900808044 ni otro SHA anterior.
