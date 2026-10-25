# Propuesta: preparar ChemUSON 0.3.0-beta.2

## Motivación

El tag anotado e inmutable `v0.3.0-beta.1` apunta al commit validado, pero el commit contiene `[skip ci]`; por ello GitHub no inició el workflow oficial tag-only. `release.yml` no admite `workflow_dispatch`, y el intento de iniciarlo manualmente fue rechazado. No hay GitHub Release ni assets oficiales para beta.1. El Build Preview de beta.1 fue temporal y no debe convertirse manualmente en publicación.

La recuperación segura es preparar la siguiente numeración monotónica, `0.3.0-beta.2`, en una rama aislada basada en `release/v0.3.0-beta.1-prep`, sin reescribir commits ni tocar el tag existente. La candidata debe pasar verificaciones focalizadas y un Build Preview nuevo. El tag beta.2 queda fuera del alcance hasta que el propietario reciba el SHA exacto y confirme el cambio de numeración.

## Alcance

- Crear `release/v0.3.0-beta.2-prep` desde el HEAD remoto verificado de `release/v0.3.0-beta.1-prep`.
- Cambiar únicamente la versión canónica, la entrada AppStream, el destino de las notas públicas que consume el workflow, la política/notas de release y los tres contratos de test que fijan la versión visible/path de notas; registrar este OpenSpec.
- Ejecutar verificaciones focalizadas, revisar que las protecciones de tags/environments y los gates de release sigan intactos, y construir un Build Preview desde el SHA candidato exacto.
- Verificar los cuatro paquetes temporales, sus hashes SHA-256 y provenance. Informar el SHA candidato y el resultado sin crear el tag beta.2.

## Fuera de alcance

- Crear, borrar, mover o recrear `v0.3.0-beta.1` o crear `v0.3.0-beta.2`.
- Crear una GitHub Release, subir assets manualmente, publicar el remoto Flatpak manualmente, cambiar permisos/gates o usar artifacts antiguos como paquetes oficiales.
- Cambiar comportamiento de producción, formato de datos, dependencias o el diseño de los workflows.
- Ejecutar la suite completa de pytest o afirmar aceptación manual específica de los nuevos artifacts sin evidencia del propietario.
