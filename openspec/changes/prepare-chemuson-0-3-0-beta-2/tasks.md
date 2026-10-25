# Tareas

- [x] Confirmar rama base, estado limpio, correcciones integradas, identidad remota del tag beta.1 y ausencia de beta.2.
- [x] Crear rama local aislada `release/v0.3.0-beta.2-prep` y registrar el baseline antes de tocar archivos de producción.
- [x] Cambiar únicamente metadatos/documentación beta.2, hacer que `release.yml` use las notas públicas beta.2 y sincronizar las tres expectativas de tests asociadas.
- [x] Ejecutar validaciones OpenSpec, compilación y pruebas focalizadas; registrar excepciones basales sin debilitarlas.
- [x] Revisar el workflow oficial, las protecciones activas de tags y los environments `beta`/`flatpak-beta`.
- [ ] Publicar la rama de preparación y ejecutar un Build Preview nuevo desde el SHA candidato exacto.
- [ ] Descargar los cuatro artifacts del Build Preview; validar identidad del SHA, provenance y SHA-256 local contra manifests.
- [ ] Dejar Git limpio e informar el SHA completo y la evidencia al propietario; no crear el tag beta.2 sin su confirmación explícita.
