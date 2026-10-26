# Diseño: beta.4 con identidad visual del proyecto y preview no publicador

## Base y aislamiento

- Base: `main`/`origin/main` en `9c171a8267d632e9bd914392680c316c28fe5e87`; `git fetch origin main` confirmó divergencia `0 0`.
- Rama: `release/v0.3.0-beta.4-prep`.
- `v0.3.0-beta.1`, `.2` y `.3` son tags anotados existentes y se preservan. No crear `v0.3.0-beta.4`.
- `gh auth status` indica que no hay sesión autenticada. La escritura remota, PR y Actions oficial quedan condicionados a credenciales accesibles; no se simulan.

## Fuente visual y variantes

- El PNG aprobado es una imagen RGB 1254×1254, PNG íntegro, SHA-256 `aab2600d2d46d87dba9dcea4266ae89ee4bbceafef0b05892b10713c310e29e6`.
- Se conserva como maestro byte a byte en `assets/branding/logo_ChemUSON.png`.
- Las variantes raster se generan por remuestreo de alta calidad sin recorte, recoloración ni cambios de proporción. Un `.ico` contiene tamaños de Windows; el SVG de distribución es un contenedor XML con el PNG original codificado en base64, no una vectorización.
- El paquete PyInstaller incluye la variante de ventana en su CArchive y usa el `.ico` como recurso PE en Windows. La aplicación carga el PNG desde la raíz de recursos de PyInstaller o desde el checkout fuente. Inno Setup utiliza el mismo `.ico` para el instalador.
- El SVG de marca sustituye el icono previo de AppImage/Flatpak, se instala bajo el ID estable existente y se declara como icono stock en AppStream. Los IDs, nombres, remotos y canales no cambian.
- README y manual presentan el emblema como identidad del proyecto ChemUSON; no como escudo oficial de Universidad de Sonora. El contraste claro/oscuro se verifica mediante rasterizado Qt de la imagen original opaca en ambos estilos y queda también sujeto a inspección visual del propietario.

## Empaquetado y limpieza

- El baseline clasifica `dist/`, `dist-appimage/`, `dist-appimage-final/`, `dist-preview/`, `dist-preview-final/` como paquetes antiguos ignorados y `build/` como intermedios PyInstaller; suman 2,499,985,408 bytes de espacio asignado. `dist-appimage-diagnostic/` contiene dos manifiestos trackeados y se conserva completo.
- No usar `git clean -fdX`: el dry-run también lista `.venv`, cachés y pycache que no se deben borrar en esta limpieza.
- Los scripts Linux existentes ya respetan `CHEMUSON_OUTPUT_DIR`. `packaging/windows/Chemuson.iss` será el único ajuste local si el test del preprocesador confirma que su `OutputDir` fijo ignora esa variable. Sin variable, se mantiene `dist-installer`, que es el contrato de los workflows Windows.
- No cambiar los comandos o rutas de los workflows oficiales por una preferencia de compilación local.

## Versión, documentación y actualizadores

- `_version.py` sigue siendo la fuente canónica. `set_version.py` actualiza versión/entrada AppStream para `0.3.0-beta.4`; las entradas históricas se conservan.
- Las notas beta.4 son bilingües y separan mantenimiento/identidad de funciones ya presentes en beta.3. Los README mantienen enlaces a la última Release que realmente existe; la vista previa no se presenta como descarga oficial.
- La API pública de GitHub confirmó en 2026-10-26 que `v0.3.0-beta.3` tiene `prerelease=false` y aparece en `/releases/latest`. Documentar como decisión excepcional del propietario. El workflow oficial todavía deriva `prerelease=true` para versiones del canal beta; no se modifica.
- El updater integrado usa SemVer de `tag_name`, no orden de fecha/título/flag, y `channel_accepts_version` excluye prereleases por versión en canal stable. La prueba focal confirma selección por canal y orden `beta.N`.
- El AppImage beta.3 publicado conserva `gh-releases-zsync|PJGV333|Chemuson|prerelease|...zsync`; el workflow de beta.4 seguirá el contrato normal `prerelease`. No cambiar preferencias del usuario. Si no hay `appimageupdatetool` local, registrar que la semántica del ejecutable AppImageUpdate no se probó dinámicamente; no atribuirle un PASS.
- Flatpak actualiza por el remoto/rama instalado (`flatpak update`), no por el updater integrado; beta.4 conserva la rama beta y el bundle preview sigue aislado.
- Proteger el downgrade y exigir checksum/provenance se mantiene sin cambios.

## Validación y publicación

- Tests focales existentes para branding, iconos, packaging, release/version y updater; compileall, Ruff focal, arquitectura, contratos estrictos OpenSpec y `git diff --check`. Colección/manifiesto node IDs sólo si la lista de tests cambia.
- Validar recursos en source y PyInstaller, AppStream, scripts, metadatos y smoke esencial de apertura/editor 2D/flujo 2D→3D en los entornos disponibles. Las pruebas visuales de Windows y temas en paquetes exactos corresponden al propietario.
- Build Preview oficial requerido: Windows portable, Windows installer, Linux AppImage Type 2 y Flatpak, mismo SHA con checksums/provenance. No llamar “preview exitoso” sin run oficial completo y revisar los cuatro artifacts.
- Detener antes del tag. No hay publicación accidental ni modificación de distribución.
