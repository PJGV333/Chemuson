# Tareas

- [x] Registrar baseline acotado y rama fuente beta.3; documentar fallos ambientales sin ejecutar pytest masivo.
- [x] Verificar merge-base, commits y diff de `origin/main`; incorporar sus README/manuales publicados a la rama beta.3 sin integrar en `main`.
- [x] Endurecer parser SemVer de prerelease/build y añadir casos inválidos relevantes.
- [x] Seleccionar releases por tag SemVer y precedencia, ignorando el flag `prerelease` como autoridad; añadir pruebas de fechas y flags contradictorios.
- [x] Añadir regresiones de beta→beta, beta→final, downgrade/equivalencia, assets incompatibles, red y checksum inválido.
- [x] Hacer dinámica la ruta versionada de release notes en el único workflow oficial y actualizar el contract estático.
- [x] Actualizar la política de versionado con rama `main`, canales legacy conservados, SemVer autoritativo, assets/notes, Latest y estrategia para beta.3.
- [x] Ejecutar únicamente pruebas focalizadas, validación workflow/OpenSpec, compileall, Ruff focal y diff-check; no ejecutar la suite completa.
- [ ] Registrar validación final, commits y SHA; publicar sólo la rama de trabajo y dejar el árbol limpio.
- [ ] Detenerse antes de merge a `main`; no crear tags, cambiar la Release beta.3 ni publicar paquetes.
