# Propuesta — revisión integral del manual de usuario contra 0.3.0-beta.4

## Problema

`docs/MANUAL_USUARIO.md` declara documentar `0.3.0-beta.4`, pero el cuerpo del manual se escribió a partir de una auditoría de octubre sobre el commit `140a080c336515650bbeea0a4e6ead67a9999b23` y se actualizó solo en la cabecera. Entre ese commit y el tag `v0.3.0-beta.4` hay 66 commits de interfaz, y el propio manual avisaba de que sus capítulos «no han sido revalidados íntegramente».

El levantamiento sistemático contra la aplicación ejecutada (`baseline.md`) encontró funciones descritas que no existen, atajos inexistentes, instrucciones falsas, duplicaciones entre capítulos y figuras que no ilustran lo que acompañan. Un manual de química con instrucciones falsas induce error científico, así que la corrección es documentalmente prioritaria.

## Alcance

- Reescribir `docs/MANUAL_USUARIO.md` describiendo **solo** funciones verificadas en el código y la ejecución de `v0.3.0-beta.4`, con clasificación explícita de lo experimental, lo pendiente y lo no disponible.
- Reorganizar la numeración y el índice: 22 capítulos en lugar de 28, sin duplicar funciones.
- Sustituir las figuras no ilustrativas por capturas auténticas de `0.3.0-beta.4` en `docs/manual-figures/`, con script reproducible y `README.md` de procedencia.
- Añadir el enlace público al manual en la sección de Ayuda y registrar como propuesta independiente la futura entrada «Ayuda → Manual de usuario».
- Actualizar `docs/manual_usuario.md` solo en su aviso de cabecera si enlaza al manual principal.

## Fuera de alcance

No se modifica **ninguna** función: ni química, ni editor molecular, ni Clean2D, ni ChemName, ni persistencia `.cmsn`, ni formatos, ni flux de publicación, ni workflows, ni instaladores, ni tags, ni releases. No se generan ni publican PDF. No se ejecuta la suite completa de pruebas: es un cambio de documentación.

Las anomalías de comportamiento detectadas (por ejemplo, `Estructura → Limpiar 2D` duplicando a `Limpiar 2D (1 paso)`, o `Archivo → Salir` sin atajo visible) **se documentan, no se arreglan**: se registran como hallazgos para que el proyecto decida.

## Módulos probables afectados

- `docs/MANUAL_USUARIO.md`
- `docs/manual-figures/` (nuevo)
- `docs/manual_usuario.md` (solo aviso de cabecera)
- `openspec/changes/revise-user-manual-beta-4/`
